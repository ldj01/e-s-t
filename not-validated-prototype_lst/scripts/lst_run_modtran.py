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

    parser.add_argument('--version',
                        action='store_true', dest='version',
                        required=False, default=False,
                        help='Reports the version of the software')

    args = parser.parse_args()

    # Report the version and exit
    if args.version:
        print(util.Version.version_text())
        sys.exit(0)  # EXIT SUCCESS

    if args.modtran_data_path is None:
        raise Exception('--modtran_data_path must be specified on the command line')

    if args.modtran_data_path == '':
        raise Exception('The MODTRAN data directory provided was empty')

    return args


TAPE5 = 'tape5'


def process_point_dir((point_path, modtran_data_path)):
    """Swim in the pool using the path

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

            finally:
                if len(output) > 0:
                    logger.info(output)
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

    if process_count > 1:
        pools = Pool(process_count)
        pools.map(process_point_dir, point_parms)
    else:
        map(process_point_dir, point_parms)


if __name__ == '__main__':
    main()
