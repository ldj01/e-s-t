#! /usr/bin/env python

'''
    FILE: extract_aux_data.py

    PURPOSE: Read the metadata XML file to determine the NARR data to be
             used.  Extract the data to specific sub-directories for
             follow-on processing.

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
from argparse import ArgumentParser
from time import sleep


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
def extract_grib_data(grb_file, pressure_numbers):
    '''
    Description:
        Configures a command line for calling the wgrib executable and then
        calls it to extract the required information from the grib file.

        The output is placed into a specified directory based on the input.
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
def retrieve_and_extract_aux_data(year, month, day, hour):
    '''
    Description:
        We are coding to use NARR data, which is provided in 3hr increments.

        Builds the strings required to retrieve the auxillary data and then
        downloads them to the current directory with specified names.

        The parameters are then extracted from the downloaded grib files.

    Note:
        We use the "Range" option in the http headers to retrieve only the
        portions of the auxillary data files that we need.
    '''

    logger = logging.getLogger(__name__)

    parms_to_extract = ['HGT', 'SPFH', 'TMP']
    # Host where the auxillary data resides
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

    for parm in parms_to_extract:
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
def process_aux_data(args):
    '''
    Description:
        Parses the XML and calls the routine to retrieve and extract the LST
        AUX data.
    '''

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

    # Retrieve the auxillary data and extract it
    retrieve_and_extract_aux_data(year, month, day, hour)


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

    if args.xml_filename == '':
        logger.fatal("No XML metadata filename provided.")
        logger.fatal("Error processing LST AUX data."
                     "  Processing will terminate.")
        sys.exit(EXIT_FAILURE)

    try:
        logger.info("Downloading and extracting LST AUX data")

        process_aux_data(args)

    except Exception, e:
        logger.exception("Error processing LST AUX data."
                         "  Processing will terminate.")
        sys.exit(EXIT_FAILURE)

    logger.info("LST AUX data downloaded and extracted")
    sys.exit(EXIT_SUCCESS)
