#! /usr/bin/env python

'''
    PURPOSE: Retieves archived NARR3D files from the CISL RDA for the dates
             requested.  Extracts the variables LST requires (HGT, TMP, SPFH)
             and repackages them into our internal location and filenames.

    PROJECT: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    LICENSE: NASA Open Source Agreement 1.3

    NOTES:

          NCEP     - National Centers for Environmental Prediction
                     http://www.ncep.noaa.gov

          NARR     - NCEP North American Regional Reanalysis

          CISL RDA - Computational & Information Systems Lab
                     Research Data Archive http://rda.ucar.edu

          NCAR     - National Center for Atmospheric Research
                     http://ncar.ucar.edu

          UCAR     - University Corporation for Atmospheric Research
                     http://www2.ucar.edu

    HISTORY:

    Date              Reason
    ----------------  --------------------------------------------------------
    July/2015         Initial implementation
'''


import os
import sys
import shutil
import logging
import errno
import commands
import requests
import calendar
import itertools
import multiprocessing as mp
from cStringIO import StringIO
from argparse import ArgumentParser
from osgeo import gdal, osr
from time import sleep
from datetime import datetime, timedelta
from contextlib import closing


from .config import settings


# ============================================================================
class Web(object):
    '''
    Description:
        Provides methods for interfacing with web resources.
    '''

    # ------------------------------------------------------------------------
    class Session(object):
        '''
        Description:
            Manages an http(s) session.
        '''

        # --------------------------------------------------------------------
        def __init__(self, max_retries=3, block_size=None, timeout=300.0):
            super(Web.Session, self).__init__()

            self.logger = logging.getLogger(__name__)

            self.session = requests.Session()

            self.timeout = timeout
            self.max_retries = max_retries

            # Determine if we are streaming or not based on block_size
            self.block_size = block_size
            self.stream = False
            if self.block_size is not None:
                self.stream = True

            adapter = requests.adapters.HTTPAdapter(max_retries=max_retries)
            self.session.mount('http://', adapter)
            self.session.mount('https://', adapter)

        # --------------------------------------------------------------------
        def login(self, login_url, login_data):
            '''
            Description:
                Provides for establishing a logged in session.
            '''

            # Login to the site
            self.session.post(url=login_url, data=login_data)

        # --------------------------------------------------------------------
        def _get_file(self, download_url, destination_file, headers=None):
            '''
            Notes: Downloading this way will place the whole source file into
                   memory before dumping to the local file.
            '''

            req = self.session.get(url=download_url, timeout=self.timeout,
                                   headers=headers)

            if not req.ok:
                self.logger.error('HTTP - Transfer of [{0}] - FAILED'
                                  .format(download_url))
                # The raise_for_status gets caught by this try's except block
                req.raise_for_status()

            # Write the downloaded data to the destination file
            with open(destination_file, 'wb') as local_fd:
                local_fd.write(req.content)

                return True

        # --------------------------------------------------------------------
        def _stream_file(self, download_url, destination_file, headers=None):
            '''
            Notes: Downloading this way streams 'block_size' of data at a
                   time.
            '''

            retrieved_bytes = 0
            with closing(self.session.get(url=download_url,
                                          timeout=self.timeout,
                                          stream=self.stream,
                                          headers=headers)) as req:
                file_size = int(req.headers['content-length'])

                # Set block size based on streaming
                if self.stream:
                    block_size = self.block_size
                else:
                    block_size = file_size

                # Write the downloaded data to the destination file
                with open(destination_file, 'wb') as local_fd:
                    for data_chunk in req.iter_content(block_size):
                        local_fd.write(data_chunk)
                        retrieved_bytes += len(data_chunk)

                if retrieved_bytes != file_size:
                    raise Exception('Transfer Failed - HTTP -'
                                    ' Retrieved {0} out of {1} bytes'
                                    .format(retrieved_bytes, file_size))
                else:
                    return True

        # --------------------------------------------------------------------
        def http_transfer_file(self, download_url, destination_file,
                               headers=None):
            '''
            Description:
                Use http to transfer a file from a source location to a
                destination file on the localhost.

            Returns:
                status_code - One of the following
                            - 200, requests.codes['ok']
                            - 404, requests.codes['not_found']:
                            - 503, requests.codes['service_unavailable']:

            Notes:
                If a 503 is returned, the logged exception should be reviewed
                to determine the real cause of the error.
            '''

            self.logger.info(download_url)

            retry_attempt = 0
            done = False
            while not done:
                status_code = requests.codes['ok']
                try:

                    done = self._stream_file(download_url,
                                             destination_file)

                    if done:
                        self.logger.info("Transfer Complete - HTTP")

                except Exception:
                    self.logger.exception('HTTP - Transfer Issue')

                    if status_code != requests.codes['not_found']:
                        if retry_attempt > self.max_retries:
                            self.logger.info('HTTP - Transfer Failed'
                                             ' - exceeded retry limit')
                            done = True
                        else:
                            retry_attempt += 1
                            sleep(int(1.5 * retry_attempt))
                    else:
                        # Not Found - So break the looping because we are done
                        done = True

            return status_code


# ============================================================================
class System(object):
    '''
    Description:
        Provides methods for interfacing with the local system.
    '''

    # ------------------------------------------------------------------------
    @staticmethod
    def execute_cmd(cmd):
        '''
        Description:
            Execute a command line and return the terminal output or raise an
            exception

        Returns:
            output - The stdout and/or stderr from the executed command.
        '''

        logger = logging.getLogger(__name__)

        output = ''

        logger.info('Executing [{0}]'.format(cmd))
        (status, output) = commands.getstatusoutput(cmd)

        if status < 0:
            message = 'Application terminated by signal [{0}]'.format(cmd)
            if len(output) > 0:
                message = ' Stdout/Stderr is: '.join([message, output])
            raise Exception(message)

        if status != 0:
            message = 'Application failed to execute [{0}]'.format(cmd)
            if len(output) > 0:
                message = ' Stdout/Stderr is: '.join([message, output])
            raise Exception(message)

        if os.WEXITSTATUS(status) != 0:
            message = ('Application [{0}] returned error code [{1}]'
                       .format(cmd, os.WEXITSTATUS(status)))
            if len(output) > 0:
                message = ' Stdout/Stderr is: '.join([message, output])
            raise Exception(message)

        return output

    # ------------------------------------------------------------------------
    @staticmethod
    def create_directory(directory):
        '''
        Description:
            Create the specified directory with some error checking.
        '''

        # Create/Make sure the directory exists
        try:
            os.makedirs(directory, mode=0755)
        except OSError as ose:
            if ose.errno == errno.EEXIST and os.path.isdir(directory):
                pass
            else:
                raise


# ============================================================================
def get_name_list(s_date, e_date):
    '''
    Description:
        Determines all of the base filenames to process based on the dates
        provided.

    Notes:
        Files typically contain 3 days.
        Special assumptions are coded for the end of the month; which
        may have 1, 2, 3, or 4 days ... depending on the month and year.
    '''

    days_3 = timedelta(days=3)
    c_date = s_date

    while c_date <= e_date:
        if c_date.day == 28:
            (x, days) = calendar.monthrange(c_date.year, c_date.month)
            yield(settings.REMOTE_NAME_FORMAT.format(c_date.year,
                                                     c_date.month,
                                                     c_date.day,
                                                     days))
            delta = timedelta(days=(days - 28 + 1))
            c_date += delta
        else:
            yield(settings.REMOTE_NAME_FORMAT.format(c_date.year,
                                                     c_date.month,
                                                     c_date.day,
                                                     c_date.day + 2))
            c_date += days_3


# ============================================================================
def process_grib_for_variable_helper(args):
    '''
    Description:
        Turns the list argument into its components for the real process
        routine, since we need multiple arguments.
    '''

    process_grib_for_variable(*args)


# ============================================================================
def process_grib_for_variable(variable, grib_file):
    '''
    Description:
        Extract the specified variable from the grib file and archive it.
    '''

    logger = logging.getLogger(__name__)

    logger.info("Processing [{0}]".format(grib_file))

    # Get the date information from the grib file
    parts = grib_file.split('.')
    year = int(parts[1][:4])
    month = int(parts[1][4:6])
    day = int(parts[1][6:8])
    hour = int(parts[1][8:])

    # Figure out the filenames to create
    hdr_name = settings.ARCHIVE_NAME_FORMAT.format(variable, year, month,
                                                   day, hour*100, 'hdr')
    grb_name = settings.ARCHIVE_NAME_FORMAT.format(variable, year, month,
                                                   day, hour*100, 'grb')

    # Create inventory/header file to extract the variable data
    cmd = ['wgrib', grib_file, '|', 'grep', variable, '>', hdr_name]
    cmd = ' '.join(cmd)
    logger.info('Executing [{0}]'.format(cmd))
    output = System.execute_cmd(cmd)
    if output is not None and len(output) > 0:
        logger.info(output)

    # Create grib files for each variable
    cmd = ['cat', hdr_name, '|',
           'wgrib', grib_file, '-i', '-grib', '-o', grb_name]
    cmd = ' '.join(cmd)
    output = ''
    logger.info('Executing [{0}]'.format(cmd))
    output = System.execute_cmd(cmd)
    if output is not None and len(output) > 0:
        logger.info(output)

    # Create new inventory/header file for the variable
    cmd = ['wgrib', grb_name, '|', 'grep', variable, '>', hdr_name]
    cmd = ' '.join(cmd)
    logger.info('Executing [{0}]'.format(cmd))
    output = System.execute_cmd(cmd)
    if output is not None and len(output) > 0:
        logger.info(output)

    # Determine the directory to place the data and create it if it does not
    # exist
    dest_path = settings.ARCHIVE_DIRECTORY_FORMAT.format(base_aux_dir, year,
                                                         month, day)
    System.create_directory(dest_path)

    # Archive the files
    logger.info('Archiving into [{0}]'.format(dest_path))
    # GRIB
    dest_file = os.path.join(dest_path, grb_name)
    shutil.copyfile(grb_name, dest_file)
    # HEADER
    dest_file = os.path.join(dest_path, hdr_name)
    shutil.copyfile(hdr_name, dest_file)

    # Cleanup the working directory
    if os.path.exists(grb_name):
        os.unlink(grb_name)
    if os.path.exists(hdr_name):
        os.unlink(hdr_name)


# ============================================================================
if __name__ == '__main__':

    # Create a command line arugment parser
    description = ('Downloads LST auxillary inputs, then archives them for'
                   ' future use.')
    parser = ArgumentParser(description=description)

    # ---- Add parameters ----
    parser.add_argument('--start-date',
                        action='store',
                        dest='start_date',
                        metavar='DATE',
                        required=False,
                        help=('The start date YYYYMMDD(inclusive)'
                              ' if requiring a range.'
                              '  Defaults to --end-date if not specified.'))

    parser.add_argument('--end-date',
                        action='store',
                        dest='end_date',
                        metavar='DATE',
                        required=False,
                        help=('The end date YYYYMMDD(inclusive)'
                              ' if requiring a range.'))

    parser.add_argument('--date',
                        action='store',
                        dest='date',
                        metavar='DATE',
                        required=False,
                        help='The date YYYYMMDD for a specific date.')

    parser.add_argument('--block-size',
                        action='store',
                        dest='block_size',
                        metavar='SIZE',
                        required=False,
                        default=16777216,  # 16MB; 33554432 = 32MB
                        help=('The block size used for streaming the download.'
                              ' Specified in bytes.'
                              '  (Default => 16777216bytes = 16MB)'))

    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s 0.0.1',
                        help='Displays the version of the software.')

    # Parse the command line parameters
    args = parser.parse_args()

    # Setup the default logger format and level. log to STDOUT
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:'
                                '%(funcName)s -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO,
                        stream=sys.stdout)

    # Get the logger
    logger = logging.getLogger(__name__)

    # Turn down the requests and urllib3 logging
    requests_logger = logging.getLogger("requests")
    requests_logger.setLevel(logging.WARNING)
    urllib3_logger = logging.getLogger("urllib3")
    urllib3_logger.setLevel(logging.WARNING)

    # Verify environment variable exists along with the directory that is
    # specified
    base_aux_dir = os.environ.get('LST_AUX_DIR')
    if base_aux_dir is None:
        logger.info('Missing environment variable LST_AUX_DIR')
        sys.exit(1)

    if not os.path.isdir(base_aux_dir):
        logger.info('LST_AUX_DIR directory does not exist')
        sys.exit(1)

    # Determine the number of threads to use, defaulting to 1
    omp_num_threads = os.environ.get('OMP_NUM_THREADS')
    if omp_num_threads is None:
        omp_num_threads = 1
    else:
        omp_num_threads = int(omp_num_threads)

    s_date = None
    e_date = None

    # Check if dates were given
    if args.date is not None:
        s_date = datetime.strptime(args.date, '%Y%m%d')
        e_date = datetime.strptime(args.date, '%Y%m%d')

    elif args.end_date is not None:
        e_date = datetime.strptime(args.end_date, '%Y%m%d')

        if args.start_date is not None:
            s_date = datetime.strptime(args.start_date, '%Y%m%d')

        else:
            s_date = e_date
    else:
        raise Exception('Must supply either --date or --end-date')

    # Start date must start on a day based on a 3day per file pattern
    day = (s_date.day - 1) / 3 * 3 + 1
    delta = timedelta(days=(day - s_date.day))
    s_date += delta

    names = list(get_name_list(s_date, e_date))

    # Establish a logged in session
    session = Web.Session(block_size=args.block_size)

    login_url = settings.ucar_login_credentials['login_url']
    login_data = settings.ucar_login_credentials['login_data']

    # Log in
    session.login(login_url, login_data)

    for name in names:
        filename = '{0}.tar'.format(name)
        logger.info('Retrieving {0}'.format(filename))

        year = name[7:11]

        url = settings.UCAR_URL_FORMAT.format(year, filename)

        session.http_transfer_file(url, filename)

        # Extract the tar'd data
        cmd = ['tar', '-xvf', filename]
        cmd = ' '.join(cmd)
        grib_files = System.execute_cmd(cmd)
        if grib_file is not None and len(grib_files) > 0:
            logger.info(grib_files)

        # Process each sub-file and archive the results using a process pool
        pool = mp.Pool(omp_num_threads)

        # For each parameter we need
        for variable in settings.NARR_VARIABLES:
            pool.map(process_grib_for_variable_helper,
                     itertools.izip(itertools.repeat(variable),
                                    grib_files.split()))

        # Cleanup - Extracted grib files
        for grib_file in grib_files.split():
            if os.path.exists(grib_file):
                os.unlink(grib_file)

        # Cleanup - The Tar ball
        if os.path.exists(filename):
            os.unlink(filename)
