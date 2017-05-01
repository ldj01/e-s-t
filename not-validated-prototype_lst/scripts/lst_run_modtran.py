#! /usr/bin/env python

'''
    File: lst_run_modtran.py

    Purpose: Runs MODTRAN on a directory structure of points with input
             information.

    Project: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    License: NASA Open Source Agreement 1.3
'''

import os
import sys
import logging
import glob
from argparse import ArgumentParser
from multiprocessing import Pool

import lst_utilities as util

from lst_grid_points import read_grid_points


def retrieve_command_line_arguments():
    """Read arguments from the command line

    Returns:
        args <arguments>: The arguments read from the command line
    """

    parser = ArgumentParser(description='Runs MODTRAN on a pre-determined'
                                        ' set of points')

    parser.add_argument('--version',
                        action='version',
                        version='Land Surface Temperature - Version 0.2.0')

    parser.add_argument('--modtran_data_path',
                        action='store', dest='modtran_data_path',
                        required=False, default=None,
                        help='Path to the MODTRAN \'DATA\' directory')

    parser.add_argument('--process_count',
                        action='store', dest='process_count',
                        required=False, default=1,
                        help='Number of processes to utilize')

    parser.add_argument('--debug',
                        action='store_true', dest='debug',
                        required=False, default=False,
                        help='Output debug messages and/or keep debug data')

    args = parser.parse_args()

    if args.modtran_data_path is None:
        raise Exception('--modtran_data_path must be specified on the'
                        'command line')

    if args.modtran_data_path == '':
        raise Exception('The MODTRAN data directory provided was empty')

    return args


class Tape6ProcessingError(Exception):
    """Exception specifically for Tape6 errors"""
    pass


TAPE6 = 'tape6'
TAPE6_SEARCH_STR = 'AREA-AVERAGED GROUND TEMPERATURE [K]'


def extract_target_surface_temp():
    """Extract area-averaged ground temperature from MODTRAN tape6 file

    Returns:
        <float>: Target surface temperature value

    Note:
        Search for the AREA-AVERAGED GROUND TEMPERATURE [K]
        Also called IMAGED-PIXEL (H2ALT) SURFACE TEMPERATURES [K]
        We are only performing the MODTRAN operation on one pixel and the
        AREA-AVERAGED is easier to extract
    """

    value = None

    with open(TAPE6, 'r') as tape6_fd:
        for line in tape6_fd:
            # Remove all whitespace and newlines
            line = ' '.join(line.strip().split())

            if line.startswith(TAPE6_SEARCH_STR):
                value = float(list(reversed(line.split()))[0])
                break

    if value is None:
        raise Tape6ProcessingError('Unable to find [{}] in tape6 file'
                                   .format(TAPE6_SEARCH_STR))

    return value


PLTOUT_ASC = 'pltout.asc'


def extract_pltout_results():
    """Extract auxiliary results
    """

    results = list()

    # Retrieve the auxillary data and extract it
    with open(PLTOUT_ASC, 'r') as pltout_fd:
        for line in pltout_fd:
            # Remove all whitespace and newlines
            line = ' '.join(line.strip().split())

            if len(line) > 0:
                results.append(line)

    return results


RESULT_HDR = 'lst_modtran.hdr'
RESULT_DATA = 'lst_modtran.data'


def create_extracted_output(tsp_value, pltout_results):
    """Creates the output file and header for the results
    """

    with open(RESULT_DATA, 'w') as radiance_data_fd:
        radiance_data_fd.write('\n'.join(pltout_results))

    with open(RESULT_HDR, 'w') as radiance_hdr_fd:
        radiance_hdr_fd.write('TARGET_PIXEL_SURFACE_TEMPERATURE {0}\n'
                              .format(tsp_value))
        radiance_hdr_fd.write('RADIANCE_RECORD_COUNT {0}\n'
                              .format(len(pltout_results)))


class ModtranProcessingError(Exception):
    """Exception specifically for MODTRAN errors"""
    pass


TAPE5 = 'tape5'


def process_point_dir((point_path, modtran_data_path)):
    """Run MODTRAN for a point and parse/format the results for later use

    Args:
        path <str>: The path to a directory containing a tape5 file
    """

    logger = logging.getLogger(__name__)

    current_directory = os.getcwd()

    # Get the real paths to simplify the logic a bit
    r_paths = [os.path.realpath(path)
               for path in glob.glob(os.path.join(point_path, '*', '*', '*'))]

    try:
        for tape5_path in r_paths:
            logger.info('Processing Directory [{}]'.format(tape5_path))
            os.chdir(tape5_path)
            # MODTRAN requires the directory to always be named 'DATA'
            util.System.create_link(modtran_data_path, 'DATA')

            output = ''
            try:
                output = util.System.execute_cmd('modtran')

                if len(output) > 0:
                    if 'STOP Error:' in output:
                        msg = ('Error processing data point [{}]'
                               .format(point_path))
                        raise ModtranProcessingError(msg)

            finally:
                if len(output) > 0:
                    logger.info(output)

            if not os.path.isfile(TAPE6):
                raise ModtranProcessingError('Missing MODTRAN output file'
                                             ' {}'.format(TAPE6))

            if not os.path.isfile(PLTOUT_ASC):
                raise ModtranProcessingError('Missing MODTRAN output file'
                                             ' {}'.format(PLTOUT_ASC))

            # Modtran is done with this point
            # So now we can parse the results and generate a specifically
            # formatted version to be used later in the processing flow
            tsp_value = extract_target_surface_temp()
            pltout_results = extract_pltout_results()

            create_extracted_output(tsp_value, pltout_results)

    finally:
        os.chdir(current_directory)


PROC_CFG_FILENAME = 'processing.conf'


def main():
    """Main processing for building the points list
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

    # Load the grid information
    (grid_points, dummy1, dummy2) = read_grid_points()

    # Cut down to just the ones we need to run MODTRAN on
    point_parms = [('{0:03}_{1:03}_{2:03}_{3:03}'.format(point.row,
                                                         point.col,
                                                         point.narr_row,
                                                         point.narr_col),
                    args.modtran_data_path)
                   for point in grid_points if point.run_modtran]

    process_count = int(args.process_count)

    try:
        if process_count > 1:
            pools = Pool(process_count)
            pools.map(process_point_dir, point_parms)
        else:
            map(process_point_dir, point_parms)
    except:
        logger.exception('Error processing points')
        raise


if __name__ == '__main__':
    main()
