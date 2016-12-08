#! /usr/bin/env python

'''
    File: lst_build_modtran_input.py

    Purpose: Builds a directory structure of points and required information
             for input to MODTRAN.

    Project: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    License: NASA Open Source Agreement 1.3
'''

import os
import sys
import logging
from argparse import ArgumentParser

from lst_grid_points import read_grid_points


class InvalidNarrDataPointError(Exception):
    """Exception for invalid NARR data points
    """
    pass


# Invalid NARR data value
INVALID_NARR_DATA_VALUE = 9.999e+20

# The number of rows and columns present in the NARR data
NARR_ROWS = 277
NARR_COLS = 349

# The number of pressure levels contained in the NARR data
PRESSURE_LAYERS = [1000, 975, 950, 925, 900,
                   875, 850, 825, 800, 775,
                   750, 725, 700, 650, 600,
                   550, 500, 450, 400, 350,
                   300, 275, 250, 225, 200,
                   175, 150, 125, 100]

'''
TODO TODO TODO TODO TODO TODO TODO TODO
TODO TODO TODO TODO TODO TODO TODO TODO
-----------------------------------------------------------------------------
I think if we know more information about the data points that will need a
NARR point, the number of ALTITUDES needed to run through MODTRAN for a point
could be reduced, further reducing the number of MODTRAN runs.

I think it would be significant for flat lands to only have to run 6 MODTRAN
runs for a point instead of 27.

Mountainous land could still end up running all 27.

I think I have the information in lst_determine_grid_points.py, and can add
it to the grid points data.  It will also need modifications with this script
and maybe the follow-on scripts and executables.
-----------------------------------------------------------------------------
TODO TODO TODO TODO TODO TODO TODO TODO
TODO TODO TODO TODO TODO TODO TODO TODO
'''
# The elevations we potentially run MODTRAN through
# We only use the ones needed based on the elevation of the NARR point
GROUND_ALT = [0.0, 0.6, 1.1,
              1.6, 2.1, 2.6,
              3.1, 3.6, 4.05]

# We can use strings for these in this code to make things simpler for
# directory creation
# The temperatures we run each elevation through
TEMPERATURES = ['273', '310', '000']
# The albedo associated with each temperature (1:1 relationship)
ALBEDOS = ['0.0', '0.0', '0.1']
# Provide them as a list of pairs
TEMP_ALBEDO_PAIRS = zip(TEMPERATURES, ALBEDOS)

# The parameters we use provided by the NARR data
# Time 0 and Time 1 Height labels and directory names
HGT_PARMS = ['HGT_t0', 'HGT_t1']
# Time 0 and Time 1 Specific Humidity labels and directory names
SPFH_PARMS = ['SPFH_t0', 'SPFH_t1']
# Time 0 and Time 1 Temperature labels and directory names
TMP_PARMS = ['TMP_t0', 'TMP_t1']
# All of them combined
PARAMETERS = HGT_PARMS + SPFH_PARMS + TMP_PARMS


def retrieve_command_line_arguments():
    """Read arguments from the command line

    Returns:
        args <arguments>: The arguments read from the command line
    """

    parser = ArgumentParser(description='Builds MODTRAN input data files for'
                                        ' a pre-determined set of points')

    parser.add_argument('--data_path',
                        action='store', dest='data_path',
                        required=False, default=None,
                        help='Specify the LST Data directory')

    parser.add_argument('--debug',
                        action='store_true', dest='debug',
                        required=False, default=False,
                        help='Output debug messages and/or keep debug data')

    parser.add_argument('--version',
                        action='store_true', dest='version',
                        required=False, default=False,
                        help='Reports the version of the software')

    args = parser.parse_args()

    # Command line arguments are required so print the help if none were
    # provided
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)  # EXIT FAILURE

    if args.version:
        sys.exit(0)  # EXIT SUCCESS

    if args.data_path is None:
        raise Exception('--data_path must be specified on the command line')

    if args.data_path == '':
        raise Exception('The LST data directory provided was empty')

    return args


def load_narr_pressure_file(parameter, layer):
    """Loads a NARR pressure layer for the parameter

    Args:
        parameter <str>: The parameter to load
        layer <str>: The pressure layer to load for the parameter

    Returns:
       [<col>[<row>]]: Two-Dimensional list of the NARR points for the
                       parameter and layer
    """

    logger = logging.getLogger(__name__)

    filename = os.path.join(parameter, '.'.join([str(layer), 'txt']))
    logger.debug('Reading NARR Pressure File [{}]'.format(filename))

    point_values = None
    with open(filename, 'r') as data_fd:
        point_values = [[float(data_fd.readline())
                         for dummy1 in xrange(NARR_ROWS)]
                        for dummy2 in xrange(NARR_COLS)]

    return point_values


def load_narr_pressure_layers(parameters, layers):
    """Loads all the parameters and presssure layes

    Args:
        parameters [<str>]: All the parameters to load
        layers [<str>]: All the pressure layers to load for each parameter

    Returns:
        <dict>: dict[parameters]->dict[layers]->[<col>[<row>]]
    """

    data = dict()
    for parameter in parameters:
        data[parameter] = dict()
        for layer in layers:
            data[parameter][layer] = dict()
            data[parameter][layer] = load_narr_pressure_file(parameter, layer)

    return data


def check_hgt(data, point, layer, time):
    """Determines"""

    # Geopotential height at time
    hgt = data[HGT_PARMS[time]][layer][point.narr_col][point.narr_row]

    if hgt == INVALID_NARR_DATA_VALUE:
        print ('HGT  : {0} {1} {2} {3} {4} {5}'
               .format(time, layer, point.narr_col, point.narr_row,
                       point.lon, point.lat))


def check_spfh(data, point, layer, time):
    """Determines"""

    # Specific humidity at time
    spfh = data[SPFH_PARMS[time]][layer][point.narr_col][point.narr_row]

    if spfh == INVALID_NARR_DATA_VALUE:
        print ('SPFH : {0} {1} {2} {3} {4} {5}'
               .format(time, layer, point.narr_col, point.narr_row,
                       point.lon, point.lat))


def check_temp(data, point, layer, time):
    """Determines"""

    # Temperature at time
    temp = data[TMP_PARMS[time]][layer][point.narr_col][point.narr_row]

    if temp == INVALID_NARR_DATA_VALUE:
        print ('TEMP : {0} {1} {2} {3} {4} {5}'
               .format(time, layer, point.narr_col, point.narr_row,
                       point.lon, point.lat))


def check_point(data, point):
    """Generate tape5 file for the current point

    Args:
        std_atmos [StdAtmosInfo]: The standard atmosphere
        data <dict>: Data structure for the parameters and pressure layers
        point <GridPointInfo>: The current point information
        interp_factor <float>: The interpolation factor to use
        head_template <str>: The template for the head of the tape5 file
        tail_template <str>: The template for the tail of the tape5 file
    """

    for layer in PRESSURE_LAYERS:
        check_hgt(data=data, point=point, layer=layer, time=0)
        check_hgt(data=data, point=point, layer=layer, time=1)
        check_spfh(data=data, point=point, layer=layer, time=0)
        check_spfh(data=data, point=point, layer=layer, time=1)
        check_temp(data=data, point=point, layer=layer, time=0)
        check_temp(data=data, point=point, layer=layer, time=1)


def determine_missing_data(grid_points):
    """
    Args:
        espa_metadata <espa.metadata>: The metadata information for the input
        data_path <str>: The directory for the NARR data files
        std_atmos [StdAtmosInfo]: The standard atmosphere
        grid_points [GridPointInfo]: List of the grid point information
    """

    # Load the NARR pressure layers into a data structure
    data = load_narr_pressure_layers(parameters=PARAMETERS,
                                     layers=PRESSURE_LAYERS)

    for point in grid_points:
        check_point(data=data, point=point)


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

    logger.info('*** Begin Determine Missing NARR Data ***')

    # Load the grid information
    (grid_points, dummy1, dummy2) = read_grid_points()

    determine_missing_data(grid_points=grid_points)

    logger.info('*** Determine Missing NARR Data - Complete ***')


if __name__ == '__main__':
    main()
