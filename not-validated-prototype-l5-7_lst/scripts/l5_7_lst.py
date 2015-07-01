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
import re
import logging
from argparse import ArgumentParser
from datetime import datetime, timedelta

# Import the metadata api found in the espa-product-formatter project
import metadata_api

# Import local modules
import lst_utilities as util

import l5_7_estimate_landsat_emissivity as estimate_landsat_emissivity
import l5_7_build_lst_data as build_lst_data


# ============================================================================
class AuxNARRGribProcessor(object):
    '''
    Description:
        Extracts parameters from the auxillary NARR data in grib format and
        places them into 'parameter' named directories.
    '''

    def __init__(self, xml_filename, base_aux_dir):
        super(AuxNARRGribProcessor, self).__init__()

        # Keep local copies of these
        self.xml_filename = xml_filename
        self.base_aux_dir = base_aux_dir

        self.parms_to_extract = ['HGT', 'SPFH', 'TMP']
        self.aux_path_template = '{0:0>4}/{1:0>2}/{2:0>2}'
        self.aux_name_template = 'narr-a_221_{0}_{1:0>2}00_000_{2}.{3}'
        self.date_template = '{0:0>4}{1:0>2}{2:0>2}'
        self.dir_template = '{0}/{1}/{2}'

        # Setup the logger to use
        self.logger = logging.getLogger(__name__)

    # ------------------------------------------------------------------------
    def extract_grib_data(self, hdr_path, grb_path, output_dir):
        '''
        Description:
            Configures a command line for calling the wgrib executable and
            then calls it to extract the required information from the grib
            file.

            The output is placed into a specified directory based on the
            input.
        '''

        util.System.create_directory(output_dir)

        with open(hdr_path, 'r') as hdr_fd:
            for line in hdr_fd.readlines():
                self.logger.debug(line.strip())
                parts = line.strip().split(':')
                record = parts[0]
                pressure = parts[6].split('=')[1]
                self.logger.debug('{0} {1}'.format(record, pressure))

                filename = '.'.join([pressure, 'txt'])
                path = os.path.join(output_dir, filename)
                cmd = ['wgrib', grb_path,
                       '-d', record,
                       '-text', '-o', path]
                cmd = ' '.join(cmd)
                self.logger.info('wgrib command = [{0}]'.format(cmd))

                # Extract the pressure data and raise any errors
                output = ''
                try:
                    output = util.System.execute_cmd(cmd)
                except Exception:
                    self.logger.error('Failed to unpack data')
                    raise
                finally:
                    if len(output) > 0:
                        self.logger.info(output)

    # ------------------------------------------------------------------------
    def extract_aux_data(self):
        '''
        Description:
            Builds the strings required to locate the auxillary data in the
            archive then extracts the parameters into parameter named
            directories.
        '''

        xml = metadata_api.parse(self.xml_filename, silence=True)
        global_metadata = xml.get_global_metadata()
        acq_date = str(global_metadata.get_acquisition_date())
        scene_center_time = str(global_metadata.get_scene_center_time())

        # Extract the individual parts from the date
        year = int(acq_date[:4])
        month = int(acq_date[5:7])
        day = int(acq_date[8:])

        # Extract the hour parts from the time and convert to an int
        hour = int(scene_center_time[:2])
        self.logger.debug('Using Acq. Date = {0} {1} {2}'
                          .format(year, month, day))
        self.logger.debug('Using Scene Center Hour = {0:0>2}'.format(hour))

        del (global_metadata)
        del (xml)

        # Determine the 3hr increments to use from the auxillary data
        # We want the one before and after the scene acquisition time
        # and convert back to formatted strings
        hour_1 = hour - (hour % 3)
        td = timedelta(hours=3)  # allows us to easily advance to the next day

        date_1 = datetime(year, month, day, hour_1)
        date_2 = date_1 + td
        self.logger.debug('Date 1 = {0}'.format(str(date_1)))
        self.logger.debug('Date 2 = {0}'.format(str(date_2)))

        for parm in self.parms_to_extract:
            # Build the source filenames for date 1
            yyyymmdd = (self.date_template
                        .format(date_1.year, date_1.month, date_1.day))
            self.logger.debug('Date 1 yyyymmdd = {0}'.format(yyyymmdd))

            hdr_1_name = (self.aux_name_template
                          .format(yyyymmdd, date_1.hour, parm, 'hdr'))
            grb_1_name = (self.aux_name_template
                          .format(yyyymmdd, date_1.hour, parm, 'grb'))
            self.logger.debug('hdr 1 = {0}'.format(hdr_1_name))
            self.logger.debug('grb 1 = {0}'.format(grb_1_name))

            tmp = (self.aux_path_template
                   .format(date_1.year, date_1.month, date_1.day))
            hdr_1_path = (self.dir_template
                          .format(base_aux_dir, tmp, hdr_1_name))
            grb_1_path = (self.dir_template
                          .format(base_aux_dir, tmp, grb_1_name))
            self.logger.info('Using {0}'.format(hdr_1_path))
            self.logger.info('Using {0}'.format(grb_1_path))

            # Build the source filenames for date 2
            yyyymmdd = (self.date_template
                        .format(date_2.year, date_2.month, date_2.day))
            self.logger.debug('Date 2 yyyymmdd = {0}'.format(yyyymmdd))

            hdr_2_name = (self.aux_name_template
                          .format(yyyymmdd, date_2.hour, parm, 'hdr'))
            grb_2_name = (self.aux_name_template
                          .format(yyyymmdd, date_2.hour, parm, 'grb'))
            self.logger.debug('hdr 2 = {0}'.format(hdr_2_name))
            self.logger.debug('grb 2 = {0}'.format(grb_2_name))

            tmp = (self.aux_path_template
                   .format(date_2.year, date_2.month, date_2.day))
            hdr_2_path = (self.dir_template
                          .format(base_aux_dir, tmp, hdr_2_name))
            grb_2_path = (self.dir_template
                          .format(base_aux_dir, tmp, grb_2_name))
            self.logger.info('Using {0}'.format(hdr_2_path))
            self.logger.info('Using {0}'.format(grb_2_path))

            # Verify that the files we need exist
            if (not os.path.exists(hdr_1_path) or
                    not os.path.exists(hdr_2_path) or
                    not os.path.exists(grb_1_path) or
                    not os.path.exists(grb_2_path)):
                raise Exception('Required LST AUX files are missing')

            output_dir = '{0}_1'.format(parm)
            self.extract_grib_data(hdr_1_path, grb_1_path, output_dir)
            output_dir = '{0}_2'.format(parm)
            self.extract_grib_data(hdr_2_path, grb_2_path, output_dir)


# ============================================================================
def generate_lst(xml_filename, base_aux_dir,
                 only_extract_aux_data=False,
                 keep_lst_temp_data=False,
                 keep_intermediate_data=False,
                 debug=False):
    '''
    Description:
        Provides the glue code for generating LST products.
    '''

    # Get the logger
    logger = logging.getLogger(__name__)

    # ------------------------------------------------------------------------
    # Retrieval and initial processing of the required AUX data
    try:
        logger.info('Extracting LST AUX data')
        current_processor = AuxNARRGribProcessor(xml_filename, base_aux_dir)
        current_processor.extract_aux_data()
    except Exception:
        logger.error('Failed processing auxillary NARR data')
        raise

    if only_extract_aux_data:
        logger.info('Stopping - User requested to stop after extracting'
                    ' LST AUX data')
        return

    # Extract the input ID from the xml filename and build some other
    # filenames
    input_id = os.path.splitext(xml_filename)[0]
    mtl_filename = '{0}_MTL.txt'.format(input_id)
    # ESPA creates the DEM for us
    dem_filename = '{0}_dem.img'.format(input_id)

    # ------------------------------------------------------------------------
    # Generate the thermal, upwelled, and downwelled radiance bands as well as
    # the atmospheric transmittance band
    cmd = ['l5_7_intermedtiate_data',
           '--xml', xml_filename,
           '--dem', dem_filename,
           '--verbose']
    if debug:
        cmd.append('--debug')

    cmd = ' '.join(cmd)
    output = ''
    try:
        logger.info('Calling [{0}]'.format(cmd))
        output = util.System.execute_cmd(cmd)
    except Exception:
        logger.error('Failed creating intermediate data')
        raise
    finally:
        if len(output) > 0:
            logger.info(output)

    # ------------------------------------------------------------------------
    # Generate Estimated Landsat Emissivity band
    try:
        current_processor = (
            estimate_landsat_emissivity.EstimateLandsatEmissivity(
                xml_filename, keep_intermediate_data))
        current_processor.generate_product()
    except Exception:
        logger.error('Failed creating Estimated Landsat Emissivity data')
        raise

    # ------------------------------------------------------------------------
    # Generate Land Surface Temperature band
    try:
        current_processor = build_lst_data.BuildLSTData(xml_filename)
        current_processor.generate_data()
    except Exception:
        logger.error('Failed processing Land Surface Temperature')
        raise

    # ------------------------------------------------------------------------
    # Cleanup
    if not keep_intermediate_data:

        # Remove the grib extraction directories
        shutil.rmtree('HGT_1', ignore_errors=True)
        shutil.rmtree('HGT_2', ignore_errors=True)
        shutil.rmtree('SPFH_1', ignore_errors=True)
        shutil.rmtree('SPFH_2', ignore_errors=True)
        shutil.rmtree('TMP_1', ignore_errors=True)
        shutil.rmtree('TMP_2', ignore_errors=True)

        # Remove the point directories generated during the core processing
        remove_dirs = set()
        point_filename = 'point_list.txt'
        with open(point_filename, 'r') as point_list_fd:
            remove_dirs = set(list([line.strip()
                                    for line in point_list_fd.readlines()]))

        for dirname in remove_dirs:
            shutil.rmtree(dirname, ignore_errors=False)

        # Finally remove the file
        os.unlink(point_filename)

    if not keep_lst_temp_data:
        util.Metadata.remove_products(xml_filename, ['lst_temp'])


# ============================================================================
if __name__ == '__main__':
    '''
    Description:
        Gathers input parameters and performs the LST processing.
    '''

    # Create a command line arugment parser
    description = ('Retrieves and generates auxillary LST inputs, then'
                   ' processes and calls other executables for LST generation')
    parser = ArgumentParser(description=description)

    # ---- Add parameters ----
    # Required parameters
    parser.add_argument('--xml',
                        action='store', dest='xml_filename',
                        required=False, default=None,
                        help='The XML metadata file to use')

    # Optional parameters
    parser.add_argument('--only-extract-aux-data',
                        action='store_true', dest='only_extract_aux_data',
                        required=False, default=False,
                        help='Stop after extracting the AUX data')

    parser.add_argument('--keep-lst-temp-data',
                        action='store_true', dest='keep_lst_temp_data',
                        required=False, default=False,
                        help=('Keep any LST temporary data added to the'
                              ' metadata XML'))

    parser.add_argument('--keep-intermediate-data',
                        action='store_true', dest='keep_intermediate_data',
                        required=False, default=False,
                        help='Keep any intermediate data generated')

    parser.add_argument('--debug',
                        action='store_true', dest='debug',
                        required=False, default=False,
                        help='Keep any debugging data')

    parser.add_argument('--version',
                        action='store_true', dest='version',
                        required=False, default=False,
                        help='Reports the version of the software')

    # Parse the command line parameters
    args = parser.parse_args()

    # Command line arguments are required so print the help if none were
    # provided
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)  # EXIT FAILURE

    # Report the version and exit
    if args.version:
        print util.Version.version_text()
        sys.exit(0)  # EXIT SUCCESS

    # Verify that the --xml parameter was specified
    if args.xml_filename is None:
        raise Exception('--xml must be specified on the command line')
        sys.exit(1)  # EXIT FAILURE

    # Setup the logging level
    log_level = logging.INFO
    if args.debug:
        log_level = logging.DEBUG

    # Setup the default logger format and level.  Log to STDOUT.
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:'
                                '%(funcName)s -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=log_level,
                        stream=sys.stdout)

    # Get the logger
    logger = logging.getLogger(__name__)

    # Verify required environment variables exists
    base_aux_dir = os.environ.get('LST_AUX_DIR')
    if base_aux_dir is None:
        logger.info('Missing environment variable LST_AUX_DIR')
        sys.exit(1)  # EXIT FAILURE

    # Not used here, only verified because a sub-executable requires it
    tmp = os.environ.get('LST_DATA_DIR')
    if tmp is None:
        logger.info('Missing environment variable LST_DATA_DIR')
        sys.exit(1)  # EXIT FAILURE

    # Not used here, only verified because a sub-executable requires it
    tmp = os.environ.get('ASTER_GED_SERVER_NAME')
    if tmp is None:
        logger.info('Missing environment variable ASTER_GED_SERVER_NAME')
        sys.exit(1)  # EXIT FAILURE

    # Verify that the base_aux_dir exists
    if not os.path.isdir(base_aux_dir):
        logger.info('LST_AUX_DIR directory does not exist')
        sys.exit(1)  # EXIT FAILURE

    # Verify that the XML filename provided is not an empty string
    if args.xml_filename == '':
        logger.fatal('No XML metadata filename provided.')
        logger.fatal('Error processing LST.  Processing will terminate.')
        sys.exit(1)  # EXIT FAILURE

    try:
        logger.info('Generating LST products')

        generate_lst(args.xml_filename, base_aux_dir,
                     args.only_extract_aux_data,
                     args.keep_lst_temp_data,
                     args.keep_intermediate_data,
                     args.debug)

    except Exception:
        logger.exception('Error processing LST.  Processing will terminate.')
        sys.exit(1)  # EXIT FAILURE

    logger.info('Completion of LST processing')
    sys.exit(0)  # EXIT SUCCESS
