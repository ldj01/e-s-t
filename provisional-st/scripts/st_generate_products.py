#! /usr/bin/env python

'''
    File: st_generate_products.py

    Purpose: Runs all of the sub-applications required to generate ST
             products.

    Project: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    License: NASA Open Source Agreement 1.3
'''

import os
import sys
import logging
import glob
import shutil
from argparse import ArgumentParser
from ConfigParser import ConfigParser

import st_utilities as util

from st_grid_points import (GRID_POINT_HEADER_NAME,
                             GRID_POINT_BINARY_NAME)

from st_build_modtran_input import PARAMETERS

import build_st_data


def retrieve_command_line_arguments():
    """Read arguments from the command line

    Returns:
        args <arguments>: The arguments read from the command line
    """

    parser = ArgumentParser(description='Creates surface temperature product')

    parser.add_argument('--version',
                        action='version',
                        version=util.Version.version_text())

    parser.add_argument('--xml',
                        action='store', dest='xml_filename',
                        required=False, default=None,
                        help='The XML metadata file to use')

    parser.add_argument('--keep-intermediate-data',
                        action='store_true', dest='intermediate',
                        required=False, default=False,
                        help='Keep any intermediate products generated')

    parser.add_argument('--keep-temporary-data',
                        action='store_true', dest='temporary',
                        required=False, default=False,
                        help='Keep any temporary files generated')

    parser.add_argument('--debug',
                        action='store_true', dest='debug',
                        required=False, default=False,
                        help='Output debug messages and/or keep debug data')

    parser.add_argument('--num_threads',
                        action='store', dest='num_threads',
                        required=False,
                        help='Number of MODTRAN processes to run at once')

    espa_group_parser = parser.add_mutually_exclusive_group(required=False)

    espa_group_parser.add_argument('--espa',
                                   action='store_true', dest='espa',
                                   required=False, default=True,
                                   help='Use ESPA processing ' +
                                   'configuration file')

    espa_group_parser.add_argument('--no_espa',
                                   action='store_false', dest='espa',
                                   required=False,
                                   help='Do not use ESPA processing ' +
                                   'configuration file')

    espa_group_parser.set_defaults(espa=True)

    args = parser.parse_args()

    # If using the espa config file, num_threads is optional
    # If not using the espa config file, provide a default num_threads
    if not args.espa and args.num_threads is None:
            args.num_threads = 1

    # Verify that the --xml parameter was specified
    if args.xml_filename is None:
        raise Exception('--xml must be specified on the command line')

    return args


def get_cfg_file_path(filename):
    """Build the full path to the config file

    Args:
        filename <str>: The name of the file to append to the full path

    Raises:
        Exception(<str>)
    """

    # Use the users home directory as the base source directory for
    # configuration
    if 'HOME' not in os.environ:
        raise Exception('[HOME] not found in environment')
    home_dir = os.environ.get('HOME')

    # Build the full path to the configuration file
    config_path = os.path.join(home_dir, '.usgs', 'espa', filename)

    return config_path


def retrieve_cfg(cfg_filename):
    """Retrieve the configuration for the cron

    Args:
        cfg_filename <str>: Name of the configuration file

    Returns:
        cfg <ConfigParser>: Configuration for ESPA cron

    Raises:
        Exception(<str>)
    """

    # Build the full path to the configuration file
    config_path = get_cfg_file_path(cfg_filename)

    if not os.path.isfile(config_path):
        raise Exception('Missing configuration file [{}]'
                        .format(config_path))

    # Create the object and load the configuration
    cfg = ConfigParser()
    cfg.read(config_path)

    return cfg


def determine_grid_points(xml_filename, data_path, debug):
    """Determines the grid points to utilize

    Args:
        xml_filename <str>: XML metadata filename
        data_path <str>: Directory for ST data files
        debug <bool>: Debug logging and processing
    """

    output = ''
    try:
        cmd = ['st_determine_grid_points.py',
               '--xml', xml_filename,
               '--data_path', data_path]

        if debug:
            cmd.append('--debug')

        output = util.System.execute_cmd(' '.join(cmd))
    finally:
        if len(output) > 0:
            logger = logging.getLogger(__name__)
            logger.info(output)


def extract_auxiliary_narr_data(xml_filename, aux_path, debug):
    """Determines the grid points to utilize

    Args:
        xml_filename <str>: XML metadata filename
        aux_path <str>: Directory for the auxiliary data files
        debug <bool>: Debug logging and processing
    """

    output = ''
    try:
        cmd = ['st_extract_auxiliary_narr_data.py',
               '--xml', xml_filename,
               '--aux_path', aux_path]

        if debug:
            cmd.append('--debug')

        output = util.System.execute_cmd(' '.join(cmd))
    finally:
        if len(output) > 0:
            logger = logging.getLogger(__name__)
            logger.info(output)


def build_modtran_input(xml_filename, data_path, debug):
    """Determines the grid points to utilize

    Args:
        xml_filename <str>: XML metadata filename
        data_path <str>: Directory for ST data files
        debug <bool>: Debug logging and processing
    """

    output = ''
    try:
        cmd = ['st_build_modtran_input.py',
               '--xml', xml_filename,
               '--data_path', data_path]

        if debug:
            cmd.append('--debug')

        output = util.System.execute_cmd(' '.join(cmd))
    finally:
        if len(output) > 0:
            logger = logging.getLogger(__name__)
            logger.info(output)


def generate_emissivity_products(xml_filename, server_name, server_path,
                                 debug):
    """Generate the required Emissivity products

    Args:
        xml_filename <str>: XML metadata filename
        server_name <str>: Name of the ASTER GED server
        server_path <str>: Path on the ASTER GED server
        debug <bool>: Debug logging and processing
    """

    output = ''
    try:
        cmd = ['estimate_landsat_emissivity.py',
               '--xml', xml_filename,
               '--aster-ged-server-name', server_name,
               '--aster-ged-server-path', server_path]

        if debug:
            cmd.append('--debug')

        output = util.System.execute_cmd(' '.join(cmd))
    finally:
        if len(output) > 0:
            logger = logging.getLogger(__name__)
            logger.info(output)

    output = ''
    try:
        cmd = ['estimate_landsat_emissivity_stdev.py',
               '--xml', xml_filename,
               '--aster-ged-server-name', server_name,
               '--aster-ged-server-path', server_path]

        if debug:
            cmd.append('--debug')

        output = util.System.execute_cmd(' '.join(cmd))
    finally:
        if len(output) > 0:
            logger = logging.getLogger(__name__)
            logger.info(output)


def run_modtran(modtran_data_path, process_count, debug):
    """Determines the grid points to utilize

    Args:
        modtran_data_path <str>: Directory for the MODTRAN 'DATA' files
        process_count <str>: Number of processes to use
        debug <bool>: Debug logging and processing
    """

    output = ''
    try:
        cmd = ['st_run_modtran.py',
               '--modtran_data_path', modtran_data_path,
               '--process_count', str(process_count)]

        if debug:
            cmd.append('--debug')

        output = util.System.execute_cmd(' '.join(cmd))
    finally:
        if len(output) > 0:
            logger = logging.getLogger(__name__)
            logger.info(output)


def generate_distance_to_cloud(xml_filename, debug):
    """Run the tool to create the distance to cloud band

    Args:
        xml_filename <str>: XML metadata filename
        debug <bool>: Debug logging and processing
    """

    output = ''
    try:
        cmd = ['st_generate_distance_to_cloud.py',
               '--xml', xml_filename]

        if debug:
            cmd.append('--debug')

        output = util.System.execute_cmd(' '.join(cmd))
    finally:
        if len(output) > 0:
            logger = logging.getLogger(__name__)
            logger.info(output)


def generate_qa(xml_filename, debug):
    """Run the tool to create the surface temperature quality band 

    Args:
        xml_filename <str>: XML metadata filename
        debug <bool>: Debug logging and processing
    """

    output = ''
    try:
        cmd = ['st_generate_qa.py',
               '--xml', xml_filename]

        if debug:
            cmd.append('--debug')

        output = util.System.execute_cmd(' '.join(cmd))
    finally:
        if len(output) > 0:
            logger = logging.getLogger(__name__)
            logger.info(output)


def convert_intermediate_bands(xml_filename, debug):
    """Run the tool to convert and scale the intermediate bands 

    Args:
        xml_filename <str>: XML metadata filename
        debug <bool>: Debug logging and processing
    """

    output = ''
    try:
        cmd = ['st_convert_bands.py',
               '--xml', xml_filename]

        if debug:
            cmd.append('--debug')

        output = util.System.execute_cmd(' '.join(cmd))
    finally:
        if len(output) > 0:
            logger = logging.getLogger(__name__)
            logger.info(output)


def cleanup_temporary_data():
    """Cleanup/remove all the ST temporary files and directories 
    """

    GRID_POINT_ELEVATION_NAME = 'grid_elevations.txt'
    MODTRAN_ELEVATION_NAME = 'modtran_elevations.txt'
    ATMOSPHERE_PARAMETERS_NAME = 'atmospheric_parameters.txt'
    USED_POINTS_NAME = 'used_points.txt'
    EMISSIVITY_HEADER_NAME = '*_emis.img.aux.xml'
    EMISSIVITY_STDEV_HEADER_NAME = '*_emis_stdev.img.aux.xml'

    # File cleanup
    cleanup_list = [GRID_POINT_HEADER_NAME, GRID_POINT_BINARY_NAME, 
                    GRID_POINT_ELEVATION_NAME, MODTRAN_ELEVATION_NAME,
                    ATMOSPHERE_PARAMETERS_NAME, USED_POINTS_NAME]

    for filename in cleanup_list:
        if os.path.exists(filename):
            os.unlink(filename)

    # Cleanup file patterns.
    cleanup_pattern_list = [EMISSIVITY_HEADER_NAME, 
                            EMISSIVITY_STDEV_HEADER_NAME]
    for pattern in cleanup_pattern_list:
        for filename in glob.glob(pattern):
            os.unlink(filename)

    # Directory cleanup
    for directory in glob.glob('[0-9][0-9][0-9]_[0-9][0-9][0-9]_'
                               '[0-9][0-9][0-9]_[0-9][0-9][0-9]'):
        shutil.rmtree(directory)

    for directory in PARAMETERS:
        if os.path.exists(directory):
            shutil.rmtree(directory)


def cleanup_intermediate_bands():
    """Cleanup/remove the intermediate bands used to make the ST band
    """

    # File cleanup
    EMISSIVITY_PATTERN = '*_emis.'
    EMISSIVITY_STDEV_PATTERN = '*_emis_stdev.'
    ATMOSPHERIC_TRANSMITTANCE_PATTERN = '*_st_atmospheric_transmittance.'
    DOWNWELLED_RADIANCE_PATTERN = '*_st_downwelled_radiance.'
    UPWELLED_RADIANCE_PATTERN = '*_st_upwelled_radiance.'
    THERMAL_RADIANCE_PATTERN = '*_st_thermal_radiance.'
    CLOUD_DISTANCE_PATTERN = '*_st_cloud_distance.'

    # Cleanup file patterns.
    cleanup_list = [EMISSIVITY_PATTERN, EMISSIVITY_STDEV_PATTERN,
                    ATMOSPHERIC_TRANSMITTANCE_PATTERN,
                    DOWNWELLED_RADIANCE_PATTERN, UPWELLED_RADIANCE_PATTERN,
                    THERMAL_RADIANCE_PATTERN, CLOUD_DISTANCE_PATTERN]

    # Only cleanup these extensions
    extensions = ('hdr', 'img')
    for extension in extensions:
        for pattern in cleanup_list:
            pattern += extension
            for filename in glob.glob(pattern):
                os.unlink(filename)


PROC_CFG_FILENAME = 'processing.conf'


def main():
    """Main processing for creating the surface temperature product 
    """

    # Command Line Arguments
    args = retrieve_command_line_arguments()

    # Check logging level
    logging_level = logging.INFO
    if args.debug:
        logging_level = logging.DEBUG

    # Setup the default logger format and level.  Log to STDOUT.
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:'
                                '%(funcName)s -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging_level,
                        stream=sys.stdout)
    logger = logging.getLogger(__name__)

    logger.info('*** Begin ST Generate Products ***')

    if (args.espa):
        # Retrieve the processing configuration
        proc_cfg = retrieve_cfg(PROC_CFG_FILENAME)

        # Determine number of processes to use
        if args.num_threads is None:
            process_count = proc_cfg.get('processing', 'omp_num_threads')
        else:
            process_count = args.num_threads

        # Determine ST data locations
        data_path = proc_cfg.get('processing', 'st_data_path')

        # Determine NARR data locations
        aux_path = proc_cfg.get('processing', 'st_aux_path')

        # Determine MODTRAN 'DATA' location
        modtran_data_path = proc_cfg.get('processing', 'modtran_data_path')

        # Determine the server name and path to get the ASTER data from
        server_name = proc_cfg.get('processing', 'aster_ged_server_name')
        server_path = proc_cfg.get('processing', 'aster_ged_server_path')
    else:
        # Determine number of processes to use
        process_count = args.num_threads

        # Determine ST data locations
        if 'ST_DATA_DIR' not in os.environ:
            raise Exception('[ST_DATA_DIR] not found in environment')
        data_path = os.environ.get('ST_DATA_DIR')

        # Determine auxiliary (NARR or MERRA2) data locations
        aux_path_variable = '{0}_AUX_DIR'.format(args.reanalysis)
        if aux_path_variable not in os.environ:
            raise Exception('[{0}] not found in environment'.
                    format(aux_path_variable))
        aux_path = os.environ.get(aux_path_variable)

        # Determine MODTRAN 'DATA' location
        if 'MODTRAN_DATA' not in os.environ:
            raise Exception('[MODTRAN_DATA] not found in environment')
        modtran_data_path = os.environ.get('MODTRAN_DATA')

        # Determine the server name and path to get the ASTER data from
        if 'ASTER_GED_SERVER' not in os.environ:
            raise Exception('[ASTER_GED_SERVER] not found in environment')
        server_name = os.environ.get('ASTER_GED_SERVER')
        if 'ASTER_GED_PATH' not in os.environ:
            raise Exception('[ASTER_GED_PATH] not found in environment')
        server_path = os.environ.get('ASTER_GED_PATH')

    # -------------- Generate the products --------------
    determine_grid_points(xml_filename=args.xml_filename,
                          data_path=data_path,
                          debug=args.debug)

    extract_auxiliary_narr_data(xml_filename=args.xml_filename,
                                aux_path=aux_path,
                                debug=args.debug)

    build_modtran_input(xml_filename=args.xml_filename,
                        data_path=data_path,
                        debug=args.debug)

    generate_emissivity_products(xml_filename=args.xml_filename,
                                 server_name=server_name,
                                 server_path=server_path,
                                 debug=args.debug)

    run_modtran(modtran_data_path=modtran_data_path,
                process_count=process_count,
                debug=args.debug)

    # Generate the thermal, upwelled, and downwelled radiance bands as well as
    # the atmospheric transmittance band
    cmd = ['st_atmospheric_parameters', '--xml', args.xml_filename]
    if args.debug:
        cmd.append('--debug')

    cmd = ' '.join(cmd)
    output = ''
    try:
        logger.info('Calling [{0}]'.format(cmd))
        output = util.System.execute_cmd(cmd)
    except Exception:
        logger.error('Failed creating atmospheric parameters and generating '
                     'intermediate data')
        raise
    finally:
        if len(output) > 0:
            logger.info(output)

    # Generate Surface Temperature band
    try:
        current_processor = build_st_data.BuildSTData(
            xml_filename=args.xml_filename)
        current_processor.generate_data()
    except Exception:
        logger.error('Failed processing Surface Temperature')
        raise

    # Build the distance to cloud band 
    generate_distance_to_cloud(xml_filename=args.xml_filename,
                               debug=args.debug)

    # Build the surface temperature quality band
    generate_qa(xml_filename=args.xml_filename,
                debug=args.debug)

    # Clean up files and directories according to user selections, or
    # for intermediate bands convert them if they are to be kept.
    if not args.temporary:
        cleanup_temporary_data()

    if args.intermediate:
        convert_intermediate_bands(xml_filename=args.xml_filename,
                                   debug=args.debug)
    else:
        cleanup_intermediate_bands()

    logger.info('*** ST Generate Products - Complete ***')


if __name__ == '__main__':
    main()
