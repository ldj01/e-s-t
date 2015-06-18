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
import re
import logging
import tempfile
import requests
from requests.exceptions import ConnectionError as Requests_ConnectionError
from requests.exceptions import Timeout as Requests_Timeout
from argparse import ArgumentParser
from time import sleep
from datetime import datetime, timedelta


# Import the metadata api found in the espa-product-formatter project
import metadata_api


# Import local modules
import lst_utilities as util


# ============================================================================
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

            util.create_directory(dir_name)

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
                    output = util.execute_cmd(cmd)
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
            byte_range  - A string containing the byte range information.
        '''

        start_bytes = list()
        end_bytes = list()
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
                        start_bytes.append(current_byte)
                        start = True

        # Combine the location values together... list of lists
        grib_ranges = list(zip(start_bytes, end_bytes))

        # Combine the ranges into a comma separated string ready for applying
        # into the Range http header parameter
        byte_range = ','.join(['{0}-{1}'.format(s, e) for s, e in grib_ranges])

        return byte_range

    # ------------------------------------------------------------------------
    def retrieve_aux_data(self, year, month, day, hour):
        '''
        Description:
            We are coding to use NARR data, which is provided in 3hr
            increments.

            Builds the strings required to retrieve the auxillary data and
            then downloads them to the current directory with specified names.

            The parameters are then extracted from the downloaded grib files.

        Returns: (dict) or None
            grb_info - A dictionary containing the filenames for each
                       parameter downloaded along with the mb numbers.
                       {
                          'PARM_NAME_1': 'NAME_OF_THE_FILE_1',
                          'PARM_NAME_2': 'NAME_OF_THE_FILE_2',
                          ...
                          'PARM_NAME_N':
                          'PARM_NAME_N': 'NAME_OF_THE_FILE_N'
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

        # Determine is the grib data exists on the server before continuing
        req = requests.head(inv_src)

        # Return to the calling code without a dict
        if not req.ok:
            return None

        # Download the inv file
        util.http_transfer_file(inv_src, inv_name)

        grb_info = dict()
        for parm in self._parms_to_extract:
            logger.info("Retrieving = {0} parameter values".format(parm))

            # Determine the specific sections of the grib file to download
            byte_range = self.determine_grib_bytes(inv_name, parm)

            # Setup the headers
            headers = {'Range': 'bytes=%s' % byte_range}

            # Figure out the local filename to create
            grb_file = grb_name.replace('.grb', '_{0}.grb'.format(parm))

            # Download the specific sections for the current parameter
            logger.info("Destination Filename = {0}".format(grb_file))
            util.http_transfer_file(grb_src, grb_file, headers=headers)

            grb_info[parm] = grb_file

        # Remove the inv file
        if os.path.exists(inv_name):
            os.unlink(inv_name)

        return grb_info

    # ------------------------------------------------------------------------
    def build_header_files(self, grb_info):
        '''
        Description:
            Configures a command line for calling the wgrib executable to
            build the headers files

        Returns: (dict)
            hdr_info - A dictionary containing the filenames for each
                       parameter downloaded along with the mb numbers.
                       {
                          'PARM_NAME_1': 'NAME_OF_THE_FILE_1',
                          'PARM_NAME_2': 'NAME_OF_THE_FILE_2',
                          ...
                          'PARM_NAME_N':
                          'PARM_NAME_N': 'NAME_OF_THE_FILE_N'
                       }
        '''

        logger = logging.getLogger(__name__)

        hdr_info = dict()
        for parm in self._parms_to_extract:

            grb_file = grb_info[parm]
            hdr_file = grb_file.replace('.grb', '.hdr')

            logger.info("Building {0} header file".format(hdr_file))

            cmd = ['wgrib', grb_file, '-h', '>', hdr_file]
            cmd = ' '.join(cmd)

            # Extract the pressure data and raise any errors
            output = ''
            try:
                logger.debug("Executing: [{0}]".format(cmd))
                output = util.execute_cmd(cmd)
            except Exception:
                logger.error("Failed reading {0} file".format(grb_file))
                raise
            finally:
                if len(output) > 0:
                    logger.info(output)

            hdr_info[parm] = hdr_file

        return hdr_info

    # ------------------------------------------------------------------------
    def archive_files(self, year, month, day, grb_info, hdr_info):
        '''
        Description:
        '''

        logger = logging.getLogger(__name__)

        if not os.path.isdir(self._base_aux_dir):
            raise Exception("Base auxillary directory does not exist")

        dest_path = '{0}/{1:0>4}/{2:0>2}/{3:0>2}'.format(self._base_aux_dir,
                                                         year, month, day)

        logger.info("Archiving into [{0}]".format(dest_path))

        util.create_directory(dest_path)

        for parm in self._parms_to_extract:

            grb_file = grb_info[parm]
            dest_file = os.path.join(dest_path, grb_file)
            shutil.copyfile(grb_file, dest_file)

            hdr_file = hdr_info[parm]
            dest_file = os.path.join(dest_path, hdr_file)
            shutil.copyfile(hdr_file, dest_file)

    # ------------------------------------------------------------------------
    def download_and_archive_aux_data(self):
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
                # 24 doesn't exist, but range doesn't include it either
                for hour in range(0, 24, 3):

                    # Retrieve the auxillary data and extract it
                    grb_info = self.retrieve_aux_data(year, month, day, hour)

                    if grb_info is None:
                        date = start_date.replace(hour=hour)
                        logger.warning("NARR data unavailable for"
                                       " {0}".format(str(date)))
                        continue

                    hdr_info = self.build_header_files(grb_info)

                    self.archive_files(year, month, day, grb_info, hdr_info)

                    # Remove the grb and hdr files
                    for parm in self._parms_to_extract:
                        if os.path.exists(grb_info[parm]):
                            os.unlink(grb_info[parm])
                        if os.path.exists(hdr_info[parm]):
                            os.unlink(hdr_info[parm])

                start_date += self._time_inc
        finally:
            # Change back to the previous directory
            os.chdir(current_directory)

            # Remove the temporary processing directory
            shutil.rmtree(temp_processing_dir, ignore_errors=True)


# ============================================================================
if __name__ == '__main__':
    '''
    Description:
        Downloads, processes, and archives NARR auxillary data for LST
        processing.
    '''

    # Create a command line arugment parser
    description = ("Downloads LST auxillary inputs, then archives them for"
                   "future use.")
    parser = ArgumentParser(description=description)

    # ---- Add parameters ----
    parser.add_argument('--start-date',
                        action='store',
                        dest='start_date',
                        required=False,
                        help=("(Optional) The starting date 'YYYYDDD'."
                              "  Required if --end-date supplied."))

    parser.add_argument('--end-date',
                        action='store',
                        dest='end_date',
                        required=False,
                        help="(Optional) The ending date 'YYYYDDD'")

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

    # Verify environment variable exists along with the directory that is
    # specified
    base_aux_dir = os.environ.get('LST_AUX_DIR')
    if base_aux_dir is None:
        logger.info("Missing environment variable LST_AUX_DIR")
        sys.exit(1)

    if not os.path.isdir(base_aux_dir):
        logger.info("LST_AUX_DIR directory does not exist")
        sys.exit(1)

    # Validate the start date is present if the end date was supplied
    if ((args.end_date is not None) and (args.start_date is None)):
        logger.error("--start-date must be specified"
                     " if --end-date is specified")
        sys.exit(1)

    # If the start date was not specified, default to the current UTC time
    if args.start_date is None:
        start_date = datetime.utcnow().replace(hour=0, minute=0,
                                               second=0, microsecond=0)
    else:
        start_date = args.start_date

    # If the end date was not specified, default to the start_date time
    if args.end_date is None:
        end_date = start_date
    else:
        end_date = args.end_date

    # If they are strings convert them to datetime objects
    if type(start_date) == str:
        length = len(start_date)
        if length < 7 or length > 10:
            logger.error("Invalid --start-date")
            sys.exit(1)

        if length == 7:
            start_date = ''.join([start_date, 'UTC'])

        try:
            start_date = datetime.strptime(start_date, '%Y%j%Z')
        except Exception:
            logger.exception("Invalid --start-date")
            sys.exit(1)

    if type(end_date) == str:
        length = len(end_date)
        if length < 7 or length > 10:
            logger.error("Invalid --end-date")
            sys.exit(1)

        if length == 7:
            end_date = ''.join([end_date, 'UTC'])

        try:
            end_date = datetime.strptime(end_date, '%Y%j%Z')
        except Exception:
            logger.exception("Invalid --end-date")
            sys.exit(1)

    if end_date < start_date:
        logger.error("--end-date must be after --start-date")
        sys.exit(1)

    try:
        logger.info("Downloading and extracting LST AUX data")

        aux_processor = AUX_Processor(base_aux_dir, start_date, end_date)

        aux_processor.download_and_archive_aux_data()

    except Exception:
        logger.exception("Error processing LST AUX data."
                         "  Processing will terminate.")
        sys.exit(1)

    logger.info("LST AUX data downloaded and extracted")
    sys.exit(0)
