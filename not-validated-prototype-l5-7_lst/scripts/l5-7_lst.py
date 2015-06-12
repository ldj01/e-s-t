#! /usr/bin/env python

'''
    FILE: l5-7_lst.py

    PURPOSE: Calls the executables required to generate the LST products.

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
import glob
import errno
import re
import commands
import logging
from argparse import ArgumentParser
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
        yyyymmdd = '{0:0>4}{1:0>2}{2:0>2}'.format(date_2.year,
                                                  date_2.month,
                                                  date_2.day)
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
        if (not os.path.exists(hdr_1_path) or
                not os.path.exists(hdr_2_path) or
                not os.path.exists(grb_1_path) or
                not os.path.exists(grb_2_path)):
            raise Exception("Required LST AUX files are missing")

        output_dir = '{0}_1'.format(parm)
        extract_grib_data(hdr_1_path, grb_1_path, output_dir)
        output_dir = '{0}_2'.format(parm)
        extract_grib_data(hdr_2_path, grb_2_path, output_dir)


# ============================================================================
def process_lst(args, base_aux_dir):
    '''
    Description:
        Provides the glue code for generating LST products.
    '''

    # Get the logger
    logger = logging.getLogger(__name__)

    # ------------------------------------------------------------------------
    # Retrieval and initial processing of the required AUX data
    try:
        logger.info("Extracting LST AUX data")

        extract_aux_data(args, base_aux_dir)
    except Exception, e:
        logger.error("Failed processing lst_download_extract_aux_data.py")
        raise e

    if args.only_extract_aux_data:
        logger.info("Stopping - User requested to stop after extracting"
                    " LST AUX data")
        return

    # Extract the product id from the xml filename and build some other
    # filenames
    product_id = os.path.splitext(args.xml_filename)[0]
    mtl_filename = '{0}_MTL.txt'.format(product_id)
    # ESPA creates the DEM for us
    dem_filename = '{0}_dem.img'.format(product_id)
    emi_filename = '{0}_emis.img'.format(product_id)

    # ------------------------------------------------------------------------
    # Generate the thermal, upwelled, and downwelled radiance bands as well as
    # the atmospheric transmittance band
    cmd = ['l5-7_intermedtiate_data',
           '--xml', args.xml_filename,
           '--dem', dem_filename,
           '--emi', emi_filename,
           '--verbose']
    if args.debug:
        cmd.append('--debug')

    cmd = ' '.join(cmd)
    output = ''
    try:
        logger.info("Calling [{0}]".format(cmd))
        output = execute_cmd(cmd)
    except Exception, e:
        logger.error("Failed processing scene_based_lst")
        raise e
    finally:
        if len(output) > 0:
            logger.info(output)

    # ------------------------------------------------------------------------
    # Generate extimated Landsat emissivity band
    cmd = ['l5-7_landsat_emissivity_from_aster_ged.py',
           '--xml', args.xml_filename]

    cmd = ' '.join(cmd)
    output = ''
    try:
        logger.info("Calling [{0}]".format(cmd))
        output = execute_cmd(cmd)
    except Exception, e:
        logger.error("Failed processing scene_based_lst")
        raise e
    finally:
        if len(output) > 0:
            logger.info(output)

    # ------------------------------------------------------------------------
    # TODO TODO TODO - Generate the LST product here
    # TODO TODO TODO - Generate the LST product here
    # TODO TODO TODO - Generate the LST product here
    # TODO TODO TODO - Generate the LST product here

    # ------------------------------------------------------------------------
    # Cleanup
    if not args.debug:

        # Remove the grib extraction directories
        shutil.rmtree('HGT_1', ignore_errors=True)
        shutil.rmtree('HGT_2', ignore_errors=True)
        shutil.rmtree('SPFH_1', ignore_errors=True)
        shutil.rmtree('SPFH_2', ignore_errors=True)
        shutil.rmtree('TMP_1', ignore_errors=True)
        shutil.rmtree('TMP_2', ignore_errors=True)

        # Remove the point directories generated during the core processing
        remove_dirs = set()
        with open('point_list.txt', "r") as point_list_fd:
            remove_dirs = set(list([line.strip()
                                    for line in point_list_fd.readlines()]))

        for dirname in remove_dirs:
            shutil.rmtree(dirname, ignore_errors=False)

        os.unlink('point_list.txt')

    # TODO TODO TODO
    # TODO TODO TODO
    # TODO TODO TODO
    # TODO TODO TODO


# ============================================================================
if __name__ == '__main__':
    '''
    Description:
        Gathers input parameters and performs the LST processing.
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
    parser.add_argument('--only-extract-aux-data',
                        action='store_true', dest='only_extract_aux_data',
                        required=False, default=False,
                        help="Stop after extracting the AUX data")

    parser.add_argument('--debug',
                        action='store_true', dest='debug',
                        required=False, default=False,
                        help="Keep any debugging data")

    # Parse the command line parameters
    args = parser.parse_args()

    # Setup the default logger format and level.  Log to STDOUT.
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:'
                                '%(funcName)s -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO,
                        stream=sys.stdout)

    # Get the logger
    logger = logging.getLogger(__name__)

    # Verify required environment variables exists
    base_aux_dir = os.environ.get('LST_AUX_DIR')
    if base_aux_dir is None:
        logger.info("Missing environment variable LST_AUX_DIR")
        sys.exit(EXIT_FAILURE)

    # Not used here, only verified because the lst executable requires it
    base_data_dir = os.environ.get('LST_DATA_DIR')
    if base_data_dir is None:
        logger.info("Missing environment variable LST_DATA_DIR")
        sys.exit(EXIT_FAILURE)

    # Verify that the base_aux_dir exists
    if not os.path.isdir(base_aux_dir):
        logger.info("LST_AUX_DIR directory does not exist")
        sys.exit(EXIT_FAILURE)

    # Verify that the XML filename provided is not an empty string
    if args.xml_filename == '':
        logger.fatal("No XML metadata filename provided.")
        logger.fatal("Error processing LST.  Processing will terminate.")
        sys.exit(EXIT_FAILURE)

    try:
        logger.info("Generating LST products")

        process_lst(args, base_aux_dir)

    except Exception, e:
        logger.exception("Error processing LST.  Processing will terminate.")
        sys.exit(EXIT_FAILURE)

    logger.info("Completion of LST processing")
    sys.exit(EXIT_SUCCESS)
