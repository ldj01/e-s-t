#! /usr/bin/env python

'''
    FILE: download_narr_aux_data.py

    PURPOSE:  Downloads NARR data.
              NCEP North American Regional Reanalysis (NARR) (32km, 25 years)

              The data is placed under the directory specified by the
              environment variable NARR_AUX_DIR.

    PROJECT: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    LICENSE: NASA Open Source Agreement 1.3

    HISTORY:

    Date              Programmer               Reason
    ----------------  ------------------------ -------------------------------
    Jan/2015          Ron Dilley               Initial implementation
'''

import os
import sys
import shutil
import errno
import re
import commands
import logging
import tempfile
import requests
from requests.exceptions import ConnectionError as Requests_ConnectionError
from requests.exceptions import Timeout as Requests_Timeout
from argparse import ArgumentParser
from time import sleep
from datetime import datetime, timedelta


# espa-common objects and methods
from espa_constants import EXIT_FAILURE
from espa_constants import EXIT_SUCCESS
import metadata_api


# ============================================================================
def execute_cmd(cmd):
    '''
    Description:
      Execute a command line and return the terminal output or raise an
      exception

    Returns:
        output - The stdout and/or stderr from the executed command.
    '''

    output = ''

    (status, output) = commands.getstatusoutput(cmd)

    if status < 0:
        message = "Application terminated by signal [%s]" % cmd
        if len(output) > 0:
            message = ' Stdout/Stderr is: '.join([message, output])
        raise Exception(message)

    if status != 0:
        message = "Application failed to execute [%s]" % cmd
        if len(output) > 0:
            message = ' Stdout/Stderr is: '.join([message, output])
        raise Exception(message)

    if os.WEXITSTATUS(status) != 0:
        message = "Application [%s] returned error code [%d]" \
                  % (cmd, os.WEXITSTATUS(status))
        if len(output) > 0:
            message = ' Stdout/Stderr is: '.join([message, output])
        raise Exception(message)

    return output


# ============================================================================
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
def http_transfer_file(download_url, destination_file, headers=None):
    '''
    Description:
        Using http transfer a file from a source location to a destination
        file on the localhost.

        HTTP headers can be specified to modify how the download happens.
    '''

    logger = logging.getLogger(__name__)

    logger.info("Transfering {0}".format(download_url))

    session = requests.Session()

    session.mount('http://', requests.adapters.HTTPAdapter(max_retries=3))
    session.mount('https://', requests.adapters.HTTPAdapter(max_retries=3))

    sleep_count = 0
    done = False
    while not done:
        if sleep_count > 3:
            raise Exception("Transfer Failed - HTTP"
                            " - exceeded retry limit")

        req = None
        try:
            req = session.get(url=download_url, headers=headers)

            if not req.ok:
                logger.error("Transfer Failed - HTTP")
                req.raise_for_status()

            with open(destination_file, 'wb') as local_fd:
                local_fd.write(req.content)

            done = True

        except:
            logger.exception("Transfer Issue - HTTP")
            sleep(int(1.5 * sleep_count))
            sleep_count = sleep_count + 1
            if sleep_count > 3:
                raise

        finally:
            if req is not None:
                req.close()

    logger.info("Transfer Complete - HTTP")


class AUX_Processor(object):
    '''
    Description:
        This provides an auxillary data processor specific to extracting
        parameters from NARR Reanalysis data for LST.

    Note:
        Without too much modification, it should be able to extract other
        parameters and also be modifiable to use MERRA as a source.
    '''

    # The base directory path to place the binary AUX files that are generated
    _base_aux_dir = None
    # The start date
    _start_date = None
    # The ending date ---- Can be equal to the start date
    _end_date = None
    # Used for advancing to the next date
    _time_inc = None

    # Host where the auxillary data resides
    _aux_hostname = None
    # Path on the host where the specifieds auxillary data is located
    _aux_path_template = None
    # Filename of the specified auxillary data on the host we require
    _aux_name_template = None

    # The parameters to extract from the grib file
    _parms_to_extract = None

    # ------------------------------------------------------------------------
    def __init__(self, base_aux_dir, start_date, end_date):
        '''
        Description:
            Initialization of the AUX_Processor object
        '''

        self._base_aux_dir = base_aux_dir
        self._start_date = start_date
        self._end_date = end_date
        self._time_inc = timedelta(days=1)

        self._aux_hostname = 'http://nomads.ncdc.noaa.gov'
        self._aux_path_template = '/data/narr/{0}/{1}/'
        self._aux_name_template = 'narr-a_221_{0}_{1:0>2}00_000.{2}'

        self._parms_to_extract = ['HGT', 'SPFH', 'TMP']

    # ------------------------------------------------------------------------
    def extract_grib_data(self, grb_info):
        '''
        Description:
            Configures a command line for calling the wgrib executable and then
            calls it to extract the required information from the grib file.

            The output is placed into a specified directory based on the input.
        '''

        logger = logging.getLogger(__name__)

        for parm in self._parms_to_extract:
            logger.info("Processing = {0} parameter values".format(parm))

            mb_numbers = grb_info[parm]['mb_numbers']
            grb_file = grb_info[parm]['filename']
            logger.debug("Grb Filename = {0}".format(grb_file))

            dir_name = \
                '{0}'.format(os.path.splitext(os.path.basename(grb_file))[0])
            logger.info("Dir Name = {0}".format(dir_name))

            create_directory(dir_name)

            index = 1
            for pressure in mb_numbers:
                filename = '.'.join([pressure, 'txt'])
                path = os.path.join(dir_name, filename)
                cmd = ['wgrib', grb_file,
                       '-d', str(index),
                       '-text',
                       '-o', path]
                cmd = ' '.join(cmd)

                # Extract the pressure data and raise any errors
                output = ''
                try:
                    logger.debug("Executing: [{0}]".format(cmd))
                    output = execute_cmd(cmd)
                except Exception, e:
                    logger.error("Failed to unpack data")
                    raise e
                finally:
                    if len(output) > 0:
                        logger.info(output)

                index = index + 1

    # ------------------------------------------------------------------------
    def determine_grib_bytes(self, inv_file, parm):
        '''
        Description:
            Reads the specified inv file and extracts the requested information
            to produce a bytes string which is compatable for using http RANGE
            for extraction of certain bytes from the source file.

        Returns: (str, list)
            bytes      - A string containing the byte range information.
            mb_numbers - A list of the pressure (mb) numbers in the order
                         found.
        '''

        start_bytes = list()
        end_bytes = list()
        mb_numbers = list()
        start = False
        with open(inv_file, 'r') as inv_fd:
            for line in inv_fd.readlines():
                line = line.rstrip('\r\n')
                split_line = line.split(':')
                current_byte = int(split_line[1])
                line_parm = split_line[3]
                if start:
                    # Subtract 1 to get the value of the previous sections
                    # ending byte
                    end_bytes.append(current_byte - 1)
                    start = False

                marker = split_line[4]
                if (line_parm == parm and re.search('\d+ mb', marker)):
                    mb_number = re.match('(\d+) mb', marker)
                    if mb_number is not None:
                        mb_numbers.append(str(mb_number.groups()[0]))
                        start_bytes.append(current_byte)
                        start = True

        # Combine the location values together... list of lists
        grib_ranges = list(zip(start_bytes, end_bytes))

        # Combine the ranges into a comma separated string ready for applying
        # into the Range http header parameter
        bytes = ','.join(['{0}-{1}'.format(s, e) for s, e in grib_ranges])

        return (bytes, mb_numbers)


    # ------------------------------------------------------------------------
    def retrieve_aux_data(self, year, month, day, hour):
        '''
        Description:
            We are coding to use NARR data, which is provided in 3hr
            increments.

            Builds the strings required to retrieve the auxillary data and
            then downloads them to the current directory with specified names.

            The parameters are then extracted from the downloaded grib files.

        Returns: (dict)
            grb_info - A dictionary containing the filenames for each
                       parameter downloaded along with the mb numbers.
                       {
                          'PARM_NAME_1':
                          {
                              'filename': 'NAME_OF_THE_FILE_1',
                              'mb_numbers': (n1, n2, ..., nn)
                          },
                          'PARM_NAME_2':
                          {
                              'filename': 'NAME_OF_THE_FILE_2',
                              'mb_numbers': (n1, n2, ..., nn)
                          },
                          ...
                          'PARM_NAME_N':
                          {
                              'filename': 'NAME_OF_THE_FILE_N',
                              'mb_numbers': (n1, n2, ..., nn)
                          }
                       }

        Note:
            We use the "Range" option in the http headers to retrieve only the
            portions of the auxillary data files that we need.
        '''

        logger = logging.getLogger(__name__)

        yyyymm = '{0}{1:0>2}'.format(year, month)
        yyyymmdd = '{0}{1:0>2}'.format(yyyymm, day)

        # Build the destination filenames
        inv_name = self._aux_name_template.format(yyyymmdd, hour, 'inv')
        grb_name = self._aux_name_template.format(yyyymmdd, hour, 'grb')

        # Build the source filenames
        data_path = self._aux_path_template.format(yyyymm, yyyymmdd)
        inv_src = '{0}{1}{2}'.format(self._aux_hostname, data_path, inv_name)
        grb_src = '{0}{1}{2}'.format(self._aux_hostname, data_path, grb_name)
        logger.debug("INV = {0}".format(inv_src))
        logger.debug("GRB = {0}".format(grb_src))

        # Download the inv file
        http_transfer_file(inv_src, inv_name)

        grb_info = dict()
        for parm in self._parms_to_extract:
            logger.info("Retrieving = {0} parameter values".format(parm))

            # Determine the specific sections of the grib file to download
            (bytes, mb_numbers) = self.determine_grib_bytes(inv_name, parm)

            headers = {'Range': 'bytes=%s' % bytes}
            grb_file = grb_name.replace('.grb', '_{0}.grb'.format(parm))

            # Download the specific sections for the current parameter
            logger.info("Destination Filename = {0}".format(grb_file))
            http_transfer_file(grb_src, grb_file, headers=headers)

            grb_info[parm] = dict()
            grb_info[parm]['filename'] = grb_file
            grb_info[parm]['mb_numbers'] = mb_numbers

        # Remove the inv file
        if os.path.exists(inv_name):
            os.unlink(inv_name)

        return grb_info

    # ------------------------------------------------------------------------
    def process_aux_data(self):
        '''
        Description:
            Parses the XML and calls the routine to retrieve and extract the
            LST AUX data.
        '''

        logger = logging.getLogger(__name__)

        temp_processing_dir = tempfile.mkdtemp(prefix='lst-aux-dir.')

        # Change to the temp directory
        current_directory = os.getcwd()
        os.chdir(temp_processing_dir)

        try:
            start_date = self._start_date
            while start_date <= self._end_date:

                # Setup the current day
                year = start_date.year
                month = start_date.month
                day = start_date.day

                # The NARR data is provided in 3hr increments
                # TODO TODO TODO - Change to the full day range
                for hour in range(0, 3, 3):

                    # Retrieve the auxillary data and extract it
                    grb_info = self.retrieve_aux_data(year, month, day, hour)

                    self.extract_grib_data(grb_info)

                start_date += self._time_inc
        finally:
            # Change back to the previous directory
            os.chdir(current_directory)

#            # Remove the temporary processing directory
#            shutil.rmtree(temp_processing_dir, ignore_errors=True)


# ============================================================================
if __name__ == '__main__':
    '''
    Description:
        Performs gathers input parameters and performs the LST processing.
    '''

    # Verify environment variable exists along with the directory that is
    # specified
    base_aux_dir = os.environ.get('NARR_AUX_DIR')
    if base_aux_dir is None:
        raise Exception("Missing environment variable NARR_AUX_DIR")

    if not os.path.isdir(base_aux_dir):
        raise Exception("NARR_AUX_DIR directory does not exist")

    # Create a command line arugment parser
    description = ("Retrieves and generates auxillary LST inputs, then"
                   " processes and calls other executables for LST generation")
    parser = ArgumentParser(description=description)

    # ---- Add parameters ----
    parser.add_argument('--start-date',
                        action='store',
                        dest='start_date',
                        required=False,
                        help="The starting date")

    parser.add_argument('--end-date',
                        action='store',
                        dest='end_date',
                        required=False,
                        help="The ending date")

    # Parse the command line parameters
    args = parser.parse_args()

    # Setup the default logger format and level. log to STDOUT.
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:'
                                '%(funcName)s -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO,
                        stream=sys.stdout)

    # Get the logger
    logger = logging.getLogger(__name__)

    # If the start date was not specified, default to the current UTC time
    if args.start_date is None:
        args.start_date = datetime.utcnow()
    # If the end date was not specified, default to the start_date time
    if args.end_date is None:
        args.end_date = args.start_date

    # If they are strings convert them to datetime objects
    if type(args.start_date) == str:
        length = len(args.start_date)
        if length < 7 or length > 10:
            logger.error("Invalid --start-date")
            sys.exit(EXIT_FAILURE)

        if length == 7:
            args.start_date = ''.join([args.start_date, 'UTC'])

        try:
            args.start_date = datetime.strptime(args.start_date, '%Y%j%Z')
        except Exception:
            logger.exception("Invalid --start-date")
            sys.exit(EXIT_FAILURE)

    if type(args.end_date) == str:
        length = len(args.end_date)
        if length < 7 or length > 10:
            logger.error("Invalid --end-date")
            sys.exit(EXIT_FAILURE)

        if length == 7:
            args.end_date = ''.join([args.end_date, 'UTC'])

        try:
            args.end_date = datetime.strptime(args.end_date, '%Y%j%Z')
        except Exception:
            logger.exception("Invalid --end-date")
            sys.exit(EXIT_FAILURE)


    try:
        logger.info("Downloading and extracting LST AUX data")

        aux_processor = AUX_Processor(base_aux_dir, args.start_date,
                                      args.end_date)

        aux_processor.process_aux_data()

    except Exception:
        logger.exception("Error processing LST AUX data."
                         "  Processing will terminate.")
        sys.exit(EXIT_FAILURE)

    logger.info("LST AUX data downloaded and extracted")
    sys.exit(EXIT_SUCCESS)
