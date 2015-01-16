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
def extract_grib_data(hdr_path, grb_path, output_dir):
    '''
    Description:
        Configures a command line for calling the wgrib executable and then
        calls it to extract the required information from the grib file.

        The output is placed into a specified directory based on the input.
    '''

    logger = logging.getLogger(__name__)

    create_directory(output_dir)

    with open(hdr_path, 'r') as hdr_fd:
        for line in hdr_fd.readlines():
            logger.debug(line.strip())
            parts = line.strip().split(':')
            record = parts[0]
            pressure = parts[6].split('=')[1]
            logger.debug("{0} {1}".format(record, pressure))

            filename = '.'.join([pressure, 'txt'])
            path = os.path.join(output_dir, filename)
            cmd = ['wgrib', grb_path,
                   '-d', record,
                   '-text', '-o', path]
            cmd = ' '.join(cmd)
            logger.info("wgrib command = [{0}]".format(cmd))

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


# ============================================================================
def extract_aux_data(args, base_aux_dir):
    '''
    Description:
        Builds the strings required to locate the auxillary data in the
        archive then extracts the parameters into paremeter named directories.
    '''

    logger = logging.getLogger(__name__)

    xml = metadata_api.parse(args.xml_filename, silence=True)
    global_metadata = xml.get_global_metadata()
    acq_date = str(global_metadata.get_acquisition_date())
    scene_center_time = str(global_metadata.get_scene_center_time())

    # Extract the individual parts from the date
    year = int(acq_date[:4])
    month = int(acq_date[5:7])
    day = int(acq_date[8:])

    # Extract the hour parts from the time and convert to an int
    hour = int(scene_center_time[:2])
    logger.debug("Using Acq. Date = {0} {1} {2}".format(year, month, day))
    logger.debug("Using Scene Center Hour = {0:0>2}".format(hour))

    del (global_metadata)
    del (xml)

    # Determine the 3hr increments to use from the auxillary data
    # We want the one before and after the scene acquisition time
    # and convert back to formatted strings
    hour_1 = hour - (hour % 3)
    td = timedelta(hours=3)  # allows us to easily advance to the next day

    date_1 = datetime(year, month, day, hour_1)
    date_2 = date_1 + td
    logger.debug("Date 1 = {0}".format(str(date_1)))
    logger.debug("Date 2 = {0}".format(str(date_2)))

    parms_to_extract = ['HGT', 'SPFH', 'TMP']
    AUX_PATH_TEMPLATE = '{0:0>4}/{1:0>2}/{2:0>2}'
    AUX_NAME_TEMPLATE = 'narr-a_221_{0}_{1:0>2}00_000_{2}.{3}'

    for parm in parms_to_extract:
        # Build the source filenames for date 1
        yyyymmdd = '{0:0>4}{1:0>2}{2:0>2}'.format(date_1.year,
                                                  date_1.month,
                                                  date_1.day)
        logger.debug("Date 1 yyyymmdd = {0}".format(yyyymmdd))

        hdr_1_name = AUX_NAME_TEMPLATE.format(yyyymmdd, date_1.hour,
                                              parm, 'hdr')
        grb_1_name = AUX_NAME_TEMPLATE.format(yyyymmdd, date_1.hour,
                                              parm, 'grb')
        logger.debug("hdr 1 = {0}".format(hdr_1_name))
        logger.debug("grb 1 = {0}".format(grb_1_name))

        tmp = AUX_PATH_TEMPLATE.format(date_1.year, date_1.month, date_1.day)
        hdr_1_path = '{0}/{1}/{2}'.format(base_aux_dir, tmp, hdr_1_name)
        grb_1_path = '{0}/{1}/{2}'.format(base_aux_dir, tmp, grb_1_name)
        logger.info("Using {0}".format(hdr_1_path))
        logger.info("Using {0}".format(grb_1_path))

        # Build the source filenames for date 2
        yyyymmdd = '{0:0>4}{1:0>2}{2:0>2}'.format(date_1.year,
                                                  date_1.month,
                                                  date_1.day)
        logger.debug("Date 2 yyyymmdd = {0}".format(yyyymmdd))

        hdr_2_name = AUX_NAME_TEMPLATE.format(yyyymmdd, date_2.hour,
                                              parm, 'hdr')
        grb_2_name = AUX_NAME_TEMPLATE.format(yyyymmdd, date_2.hour,
                                              parm, 'grb')
        logger.debug("hdr 2 = {0}".format(hdr_2_name))
        logger.debug("grb 2 = {0}".format(grb_2_name))

        tmp = AUX_PATH_TEMPLATE.format(date_2.year, date_2.month, date_2.day)
        hdr_2_path = '{0}/{1}/{2}'.format(base_aux_dir, tmp, hdr_2_name)
        grb_2_path = '{0}/{1}/{2}'.format(base_aux_dir, tmp, grb_2_name)
        logger.info("Using {0}".format(hdr_2_path))
        logger.info("Using {0}".format(grb_2_path))

        # Verify that the files we need exist
        if (not os.path.exists(hdr_1_path)
                or not os.path.exists(hdr_2_path)
                or not os.path.exists(grb_1_path)
                or not os.path.exists(grb_2_path)):
            raise Exception("Required LST AUX files are missing")

        output_dir = '{0}_1'.format(parm)
        extract_grib_data(hdr_1_path, grb_1_path, output_dir)
        output_dir = '{0}_2'.format(parm)
        extract_grib_data(hdr_2_path, grb_2_path, output_dir)


# ============================================================================
if __name__ == '__main__':
    '''
    Description:
        Performs gathers input parameters and performs the LST processing.
    '''

    # Create a command line arugment parser
    description = ("Extracts auxillary LST inputs from the auxillary archive"
                   "for the specified data.")
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

    # Verify environment variable exists along with the directory that is
    # specified
    base_aux_dir = os.environ.get('LST_AUX_DIR')
    if base_aux_dir is None:
        logger.info("Missing environment variable LST_AUX_DIR")
        sys.exit(EXIT_FAILURE)

    if not os.path.isdir(base_aux_dir):
        logger.info("LST_AUX_DIR directory does not exist")
        sys.exit(EXIT_FAILURE)

    if args.xml_filename == '':
        logger.fatal("No XML metadata filename provided.")
        logger.fatal("Error processing LST AUX data."
                     "  Processing will terminate.")
        sys.exit(EXIT_FAILURE)

    try:
        logger.info("Extracting LST AUX data")

        extract_aux_data(args, base_aux_dir)

    except Exception, e:
        logger.exception("Error processing LST AUX data."
                         "  Processing will terminate.")
        sys.exit(EXIT_FAILURE)

    logger.info("LST AUX data extracted")
    sys.exit(EXIT_SUCCESS)
