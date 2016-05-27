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
import math
import array
import struct
import numpy as np
from argparse import ArgumentParser
from collections import namedtuple
from osgeo import gdal, osr

from espa import Metadata
from lst_exceptions import MissingBandError
import lst_utilities as util

from lst_grid_points import PointInfo, write_grid_points, read_grid_points


def retrieve_command_line_arguments():
    """Read arguments from the command line

    Returns:
        args <arguments>: The arguments read from the command line
    """

    # Create a command line arugement parser
    parser = ArgumentParser(description='Runs MODTRAN on a pre-determined'
                                        ' set of points')

    # ---- Add Arguments ----
    parser.add_argument('--xml',
                        action='store', dest='xml_filename',
                        required=False, default=None,
                        help='The XML metadata file to use')

    parser.add_argument('--debug',
                        action='store_true', dest='debug',
                        required=False, default=False,
                        help='Output debug messages and/or keep debug data')

    parser.add_argument('--version',
                        action='store_true', dest='version',
                        required=False, default=False,
                        help='Reports the version of the software')

    # Parse the command line arguments
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
        raise Exception('The XML metadata filename provided was empty')

    return args


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

    # XML Metadata
    espa_metadata = Metadata()
    espa_metadata.parse(xml_filename=args.xml_filename)

    raise NotImplementedError('I\'m not implemented yet!!!!!!!!!!!!')


if __name__ == '__main__':
    main()
