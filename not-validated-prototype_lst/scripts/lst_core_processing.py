#! /usr/bin/env python

'''
    FILE: lst_core_processing.py

    PURPOSE: Calls the executables required to generate the LST products.

    PROJECT: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    LICENSE: NASA Open Source Agreement 1.3
'''

import os
import sys
import shutil
import logging
from argparse import ArgumentParser

# Import local modules
import lst_utilities as util
from lst_environment import Environment

from extract_auxiliary_narr_data import AuxNARRGribProcessor
import estimate_landsat_emissivity
import build_lst_data


def generate_lst(xml_filename,
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

    # Retrieval and initial processing of the required AUX data
    try:
        logger.info('Extracting LST AUX data')
        current_processor = AuxNARRGribProcessor(xml_filename)
        current_processor.extract_aux_data()
    except Exception:
        logger.error('Failed processing auxillary NARR data')
        raise

    if only_extract_aux_data:
        logger.info('Stopping - User requested to stop after extracting'
                    ' LST AUX data')
        return

    # Generate the thermal, upwelled, and downwelled radiance bands as well as
    # the atmospheric transmittance band
    cmd = ['lst_intermediate_data',
           '--xml', xml_filename,
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

    # Generate Estimated Landsat Emissivity band
    try:
        current_processor = (
            estimate_landsat_emissivity.EstimateLandsatEmissivity(
                xml_filename, keep_intermediate_data))
        current_processor.generate_product()
    except Exception:
        logger.error('Failed creating Estimated Landsat Emissivity data')
        raise

    # Generate Land Surface Temperature band
    try:
        current_processor = build_lst_data.BuildLSTData(xml_filename)
        current_processor.generate_data()
    except Exception:
        logger.error('Failed processing Land Surface Temperature')
        raise

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
        print(util.Version.version_text())
        sys.exit(0)  # EXIT SUCCESS

    # Verify that the --xml parameter was specified
    if args.xml_filename is None:
        raise Exception('--xml must be specified on the command line')

    # Verify that the XML filename provided is not an empty string
    if args.xml_filename == '':
        raise Exception('No XML metadata filename provided.')

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

    try:
        environment = Environment()
        del environment

        logger.info('Generating LST products')
        generate_lst(args.xml_filename,
                     args.only_extract_aux_data,
                     args.keep_lst_temp_data,
                     args.keep_intermediate_data,
                     args.debug)

    except Exception:
        logger.exception('Error processing LST.  Processing will terminate.')
        sys.exit(1)  # EXIT FAILURE

    logger.info('Completion of LST processing')
    sys.exit(0)  # EXIT SUCCESS
