#! /usr/bin/env python

'''
    FILE: do_lst.py

    PURPOSE: Performs retrieval and setup of required auxillary inputs.
             Calls the MODTRAN and other C-executables required to generate
             the LST products.

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
import errno
import re
import commands
import logging
import requests
from argparse import ArgumentParser

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
    '''

    logger = logging.getLogger(__name__)

    logger.info("Transfering {0}".format(download_url))

    req = requests.get(download_url, headers=headers)

    if not req.ok:
        logger.error("Transfer Failed - HTTP")
        req.raise_for_status()

    try:
        with open(destination_file, 'wb') as local_fd:
            local_fd.write(req.content)
    except:
        logger.error("Transfer Failed - HTTP")
        raise
    finally:
        req.close()
    logger.info("Transfer Complete - HTTP")


# ============================================================================
def determine_grib_bytes(inv_file, parm):
    '''
    Description:
        Reads the specified inv file and extracts the requested information
        to product a bytes string which is compatable for using http RANGE
        for extraction of certain bytes from the source file.

    Returns: (str, list)
        bytes - A string containing the byte range information.
        pressure_numbers - A list of the pressure numbers in the order found.
    '''

    start_bytes = list()
    end_bytes = list()
    pressure_numbers = list()
    start = False
    with open(inv_file, 'r') as inv_fd:
        for line in inv_fd.readlines():
            line = line.rstrip('\r\n')
            split_line = line.split(':')
            current_byte = int(split_line[1])
            line_parm = split_line[3]
            if start:
                end_bytes.append(current_byte - 1)
                start = False

            marker = split_line[4]
            if (line_parm == parm and re.search('\d+ mb', marker)):
                pressure_number = re.match('(\d+) mb', marker)
                if pressure_number is not None:
                    pressure_numbers.append(str(pressure_number.groups()[0]))
                    start_bytes.append(current_byte)
                    start = True

    grib_ranges = list(zip(start_bytes, end_bytes))

    bytes = ','.join(['{0}-{1}'.format(s, e) for s, e in grib_ranges])
    return (bytes, pressure_numbers)


# ============================================================================
def extract_grib_data(grb_file, pressure_numbers):
    '''
    Description:
        TODO TODO TODO
    '''

    logger = logging.getLogger(__name__)

    dir_name = '{0}'.format(os.path.splitext(os.path.basename(grb_file))[0])

    create_directory(dir_name)

    index = 1
    for pressure in pressure_numbers:
        filename = '.'.join([pressure, 'txt'])
        path = os.path.join(dir_name, filename)
        cmd = ['wgrib', grb_file, '-d', str(index), '-text', '-o', path]
        cmd = ' '.join(cmd)

        # Extract the pressure data and raise any errors
        output = ''
        try:
            output = execute_cmd(cmd)
        except Exception, e:
            logger.error("Failed to unpack data")
            raise e
        finally:
            if len(output) > 0:
                logger.info(output)

        index = index + 1


# ============================================================================
def retrieve_aux_data(year, month, day, hour):
    '''
    Description:
        We are coding to use NARR data, which is provided in 3hr increments.

        Builds the strings required to retrieve the auxillary data and then
        downloads them to the current directory with specified names.

    Note:
        We use the "Range" option in the http headers to retrieve only the
        portions of the auxillary data files that we need.
    '''

    logger = logging.getLogger(__name__)

    valid_parms = ['HGT', 'SPFH', 'TMP']
    # Host where the auxillary data resides
# TODO TODO TODO - Maybe make this an environment variable
    AUX_HOSTNAME = 'http://nomads.ncdc.noaa.gov'
    # Path on the host where the specifieds auxillary data is located
    AUX_PATH_TEMPLATE = '/data/narr/{0}/{1}/'
    # Filename of the specified auxillary data we require
    AUX_NAME_TEMPLATE = 'narr-a_221_{0}_{1}00_000.{2}'

    # Determine the 3hr increments to use from the auxillary data
    # We want the one before and after the scene acquisition time
    # and convert back to formatted strings
    hour_1 = '{0:0>2}'.format(hour - (hour % 3))
    hour_2 = '{0:0>2}'.format(hour + (3-(hour % 3)))

    logger.debug("HOUR 1 = {0}".format(hour_1))
    logger.debug("HOUR 2 = {0}".format(hour_2))

    yyyymm = '{0}{1}'.format(year, month)
    yyyymmdd = '{0}{1}'.format(yyyymm, day)

    # Build the destination filenames
    inv_1_name = AUX_NAME_TEMPLATE.format(yyyymmdd, hour_1, 'inv')
    grb_1_name = AUX_NAME_TEMPLATE.format(yyyymmdd, hour_1, 'grb')
    inv_2_name = AUX_NAME_TEMPLATE.format(yyyymmdd, hour_2, 'inv')
    grb_2_name = AUX_NAME_TEMPLATE.format(yyyymmdd, hour_2, 'grb')

    # Build the source filenames
    data_path = AUX_PATH_TEMPLATE.format(yyyymm, yyyymmdd)
    inv_1_src = '{0}{1}{2}'.format(AUX_HOSTNAME, data_path, inv_1_name)
    grb_1_src = '{0}{1}{2}'.format(AUX_HOSTNAME, data_path, grb_1_name)
    inv_2_src = '{0}{1}{2}'.format(AUX_HOSTNAME, data_path, inv_2_name)
    grb_2_src = '{0}{1}{2}'.format(AUX_HOSTNAME, data_path, grb_2_name)
    logger.debug("INV 1 = {0}".format(inv_1_src))
    logger.debug("GRB 1 = {0}".format(grb_1_src))
    logger.debug("INV 2 = {0}".format(inv_2_src))
    logger.debug("GRB 2 = {0}".format(grb_2_src))

    # Download the inv files
    http_transfer_file(inv_1_src, inv_1_name)
    http_transfer_file(inv_2_src, inv_2_name)

    for parm in valid_parms:
        logger.info("Retrieving = {0} parameters for time 1".format(parm))
        # Determine the specific sections of the grib file to download
        (bytes, pressure_numbers) = determine_grib_bytes(inv_1_name, parm)
        headers = {'Range': 'bytes=%s' % bytes}
        grb_file = '{0}_1.grb'.format(parm)
        # Download the specific sections
        http_transfer_file(grb_1_src, grb_file, headers=headers)
        # Extract the sections to text files
        extract_grib_data(grb_file, pressure_numbers)
        if os.path.exists(grb_file):
            os.unlink(grb_file)

        logger.info("Retrieving = {0} parameters for time 2".format(parm))
        # Determine the specific sections of the grib file to download
        (bytes, pressure_numbers) = determine_grib_bytes(inv_2_name, parm)
        headers = {'Range': 'bytes=%s' % bytes}
        grb_file = '{0}_2.grb'.format(parm)
        # Download the specific sections
        http_transfer_file(grb_2_src, grb_file, headers=headers)
        # Extract the sections to text files
        extract_grib_data(grb_file, pressure_numbers)
        if os.path.exists(grb_file):
            os.unlink(grb_file)

    # Remove the inv files
    if os.path.exists(inv_1_name):
        os.unlink(inv_1_name)
    if os.path.exists(inv_2_name):
        os.unlink(inv_2_name)


# ============================================================================
def process_lst(args):
    '''
    Description:
        TODO TODO TODO
    '''

    # get the logger
    logger = logging.getLogger(__name__)

    xml = metadata_api.parse(args.xml_filename, silence=True)
    global_metadata = xml.get_global_metadata()
    acq_date = str(global_metadata.get_acquisition_date())
    scene_center_time = str(global_metadata.get_scene_center_time())

    # Extract the individual parts from the date
    year = acq_date[:4]
    month = acq_date[5:7]
    day = acq_date[8:]

    # Extract the hour parts from the time and convert to an int
    hour = int(scene_center_time[:2])
    logger.debug("Using Acq. Date = {0} {1} {2}".format(year, month, day))
    logger.debug("Using Scene Center Hour = {0:0>2}".format(hour))

    del (global_metadata)
    del (xml)

    # Retrieve the auxillary data
    retrieve_aux_data(year, month, day, hour)

    # Check for stop flag
    if args.stop_at_narr:
        logger.info("User requested to stop after narr processing.")
        return

    # Tape 5 generation
    # TODO TODO TODO

    # Check for stop flag
    if args.stop_at_tape5:
        logger.info("User requested to stop after tape5 generation.")
        return


    # TODO TODO TODO
    # TODO TODO TODO
    # TODO TODO TODO
    # TODO TODO TODO
    # TODO TODO TODO


# ============================================================================
if __name__ == '__main__':
    '''
    Description:
        Performs gathers input parameters and performs the LST processing.
    '''

    # Create a command line arugment parser
    description = ("Retrieves and generates auxillary LST inputs, then"
                   " processes and calls other executables for LST generation")
    parser = ArgumentParser(description=description)

    # ---- Add parameters ----
    # Required parameters
    parser.add_argument('--xml',
                        action='store', dest='xml_filename', required=True,
                        help="The XML metadata file to use")

    # Optional parameters
    parser.add_argument('--stop-at-narr',
                        action='store_true', dest='stop_at_narr',
                        required=False, default=False,
                        help="Stop after NARR retrieval and processing")

    parser.add_argument('--stop-at-tape5',
                        action='store_true', dest='stop_at_tape5',
                        required=False, default=False,
                        help="Stop after tape5 generation")

    # Parse the command line parameters
    args = parser.parse_args()

    # setup the default logger format and level. log to STDOUT.
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:'
                                '%(funcName)s -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO)

    # get the logger
    logger = logging.getLogger(__name__)

    logger.info("Generating LST products")

    if args.xml_filename == '':
        logger.fatal("No XML metadata filename provided.")
        logger.fatal("Error processing LST.  Processing will terminate.")
        sys.exit(EXIT_FAILURE)

    try:
        process_lst(args)

    except Exception, e:
        logger.fatal(str(e))
        logger.fatal("Error processing LST.  Processing will terminate.")
        sys.exit(EXIT_FAILURE)

    logger.info("Completion of LST")
    sys.exit(EXIT_SUCCESS)
