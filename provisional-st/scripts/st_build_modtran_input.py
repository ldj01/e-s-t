#! /usr/bin/env python3

'''
    File: st_build_modtran_input.py

    Purpose: Builds a directory structure of points and required information
             for input to MODTRAN.

    Project: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    License: NASA Open Source Agreement 1.3
'''

import os
import sys
import logging
import math
from argparse import ArgumentParser
from collections import namedtuple
from bisect import bisect
import numpy as np

from espa import Metadata
import st_utilities as util

from st_grid_points import read_grid_points

GRID_ELEVATION_NAME = 'grid_elevations.txt'
MODTRAN_ELEVATION_NAME = 'modtran_elevations.txt'


class InvalidDataPointError(Exception):
    """Exception for invalid data points
    """
    pass


# Invalid data value.  This is used in NARR.  It is a different value in 
# MERRA and GEOS5, but when we write the MERRA  or GEOS5 NetCDF values to the
# ASCII files, instead of "--" we write this value for invalid data.
INVALID_DATA_VALUE = 9.999e+20

# The number of rows and columns present in the NARR data
NARR_ROWS = 277
NARR_COLS = 349
# The number of rows and columns present in the MERRA2 and GEOS5 data
MERRA_ROWS = 361 
MERRA_COLS = 576

# The pressure levels contained in the NARR data
NARR_PRESSURE_LAYERS = [1000, 975, 950, 925, 900,
                        875, 850, 825, 800, 775,
                        750, 725, 700, 650, 600,
                        550, 500, 450, 400, 350,
                        300, 275, 250, 225, 200,
                        175, 150, 125, 100]

# The pressure levels contained in the MERRA2 and GEOS5 data
MERRA_PRESSURE_LAYERS = [1000, 975, 950, 925, 900, 875, 850, 825, 800,
                         775, 750, 725, 700, 650, 600, 550, 500, 450,
                         400, 350, 300, 250, 200, 150, 100, 70, 50, 40,
                         30, 20, 10, 7, 5, 4, 3, 2, 1, 0.7, 0.5, 0.4, 0.3, 0.1]

# The elevations we potentially run MODTRAN through
# We only use the ones needed based on the elevation of the reanalysis point
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
TEMP_ALBEDO_PAIRS = list(zip(TEMPERATURES, ALBEDOS))

# The parameters we use provided by the reanalysis data
# Time 0 and Time 1 Height labels and directory names
NARR_HGT_PARMS = ['HGT_t0', 'HGT_t1']
MERRA_HGT_PARMS = ['H_t0', 'H_t1']
# Time 0 and Time 1 Specific Humidity labels and directory names
NARR_SPFH_PARMS = ['SPFH_t0', 'SPFH_t1']
MERRA_SPFH_PARMS = ['QV_t0', 'QV_t1']
# Time 0 and Time 1 Temperature labels and directory names
NARR_TMP_PARMS = ['TMP_t0', 'TMP_t1']
MERRA_TMP_PARMS = ['T_t0', 'T_t1']
# All of them combined
NARR_PARAMETERS = NARR_HGT_PARMS + NARR_SPFH_PARMS + NARR_TMP_PARMS
MERRA_PARAMETERS = MERRA_HGT_PARMS + MERRA_SPFH_PARMS + MERRA_TMP_PARMS

# MODTRAN template files located in the static data directory
MODTRAN_HEAD = 'modtran_head.txt'
MODTRAN_TAIL = 'modtran_tail.txt'
# Name of the file to write for MODTRAN to execute on
TAPE5 = 'tape5'


def retrieve_command_line_arguments():
    """Read arguments from the command line

    Returns:
        args <arguments>: The arguments read from the command line
    """

    parser = ArgumentParser(description='Builds MODTRAN input data files for'
                                        ' a pre-determined set of points')

    parser.add_argument('--version',
                        action='version',
                        version=util.Version.version_text())

    parser.add_argument('--xml',
                        action='store', dest='xml_filename',
                        required=False, default=None,
                        help='The XML metadata file to use')

    parser.add_argument('--data_path',
                        action='store', dest='data_path',
                        required=False, default=None,
                        help='Specify the ST Data directory')

    parser.add_argument('--reanalysis',
                        action='store', dest='reanalysis',
                        required=False, default='MERRA2',
                        choices=['NARR','MERRA2','GEOS5'],
                        help='Reanalysis source - NARR, MERRA2 or GEOS5')

    parser.add_argument('--debug',
                        action='store_true', dest='debug',
                        required=False, default=False,
                        help='Output debug messages and/or keep debug data')

    args = parser.parse_args()

    # Command line arguments are required so print the help if none were
    # provided
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)  # EXIT FAILURE

    if args.xml_filename is None:
        raise Exception('--xml must be specified on the command line')

    if args.xml_filename == '':
        raise Exception('The XML metadata filename provided was empty')

    if args.data_path is None:
        raise Exception('--data_path must be specified on the command line')

    if args.data_path == '':
        raise Exception('The ST data directory provided was empty')

    return args


def load_modtran_template_file(data_path, filename):
    """Loads the specified MODTRAN template file

    Args:
        data_path <str>: The directory for ST data files
        filename <str>: The filename to load

    Args:
        <str>: The contents of the requested file
    """

    filename = os.path.join(data_path, filename)
    file_data = None

    with open(filename, 'r') as file_fd:
        file_data = file_fd.read()

    return file_data


def load_pressure_file(parameter, layer, reanalysis):
    """Loads a reanalysis pressure layer for the parameter

    Args:
        parameter <str>: The parameter to load
        layer <str>: The pressure layer to load for the parameter
        reanalysis <str>: Reanalysis source: NARR, MERRA2 or GEOS5

    Returns:
       [<col>[<row>]]: Two-Dimensional list of the reanalysis points for the
                       parameter and layer
    """

    logger = logging.getLogger(__name__)

    if reanalysis == "NARR":
        file_ext = 'txt'
    else:
        file_ext = 'npy'
    filename = os.path.join(parameter, '.'.join([str(layer), file_ext]))
    logger.debug('Reading Pressure File [{}]'.format(filename))

    if reanalysis == "NARR":
        point_values = None
        with open(filename, 'r') as data_fd:
            point_values = [[float(data_fd.readline())
                             for dummy1 in range(NARR_COLS)]
                            for dummy2 in range(NARR_ROWS)]
    else:
        point_values = np.load(filename)

    return point_values


def load_pressure_layers(parameters, layers, reanalysis):
    """Loads all the parameters and presssure layers

    Args:
        parameters [<str>]: All the parameters to load
        layers [<str>]: All the pressure layers to load for each parameter
        reanalysis <str>: Reanalysis source: NARR, MERRA2 or GEOS5

    Returns:
        <dict>: dict[parameters]->dict[layers]->[<col>[<row>]]
    """

    data = dict()
    for parameter in parameters:
        data[parameter] = dict()
        for layer in layers:
            data[parameter][layer] = dict()
            data[parameter][layer] = load_pressure_file(parameter, layer,
                                                        reanalysis)

    return data


def determine_interp_factor(value, before, after):
    """Determine interpolation factor for removal of division in follow-on
       code

    For the formula:
        v = v0 + (v1 - v0) * ((t - t0) / (t1 - t0))

    I'm defining the interpolation factor as this component of the formula:

        ((t - t0) / (t1 - t0))

    Args:
        value <number>: The value to interpolate for
        before <number>: The value before
        after <number>: The value after

    Returns:
        <float>: The interpolation factor
    """

    f_value = float(value)
    f_before = float(before)
    f_after = float(after)

    return (f_value - f_before) / (f_after - f_before)


STD_ATMOS_FILENAME = 'std_mid_lat_summer_atmos.txt'
StdAtmosInfo = namedtuple('StdAtmosInfo',
                          ('hgt', 'pressure', 'temp', 'rh'))


def load_std_atmosphere(data_path):
    """Loads the standard atmosphere into an iterable

    Args:
        data_path <str>: The directory for ST data files

    Returns:
        <iterable>: Essentially a list of StdAtmosInfo objects
    """

    logger = logging.getLogger(__name__)

    filename = os.path.join(data_path, STD_ATMOS_FILENAME)
    logger.debug('Reading Standard Atmosphere File [{}]'.format(filename))

    with open(filename, 'r') as data_fd:
        for line in data_fd.readlines():
            (hgt, pressure, temp, rel_hj) = line.strip().split()
            yield(StdAtmosInfo(hgt=float(hgt), pressure=float(pressure),
                               temp=float(temp), rh=float(rel_hj)))


EQUATORIAL_RADIUS = 6378137.0
POLAR_RADIUS = 6356752.3142
STANRDARD_GRAVITY_IN_M_PER_SEC_SQRD = 9.80665
EQUATORIAL_RADIUS_IN_KM = EQUATORIAL_RADIUS / 1000.0
POLAR_RADIUS_IN_KM = POLAR_RADIUS / 1000.0
INV_R_MIN_SQRD = 1.0 / (POLAR_RADIUS_IN_KM * POLAR_RADIUS_IN_KM)
INV_R_MAX_SQRD = 1.0 / (EQUATORIAL_RADIUS_IN_KM * EQUATORIAL_RADIUS_IN_KM)
INV_STD_GRAVITY = 1.0 / STANRDARD_GRAVITY_IN_M_PER_SEC_SQRD


def determine_geometric_hgt(data, point, layer, time, reanalysis):
    """Determines geometric height for a layer in the point at time x

    Geopotential height is converted to geometric height for both t0 and t1
    heights

    Args:
        data <dict>: Data structure for the parameters and pressure layers
        point <GridPointInfo>: The current point information
        layer <str>: The layer to get the height for
        time <int>: Time 0 or Time 1
        reanalysis <str>: Reanalysis source: NARR, MERRA2 or GEOS5

    Returns:
        <float>: The geometric height for the layer
    """

    # Geopotential height at time
    if reanalysis == "NARR":
        hgt = data[NARR_HGT_PARMS[time]][layer]\
                [point.reanalysis_row-1][point.reanalysis_col-1]
        # For NARR, invalid data is an error case.
        if abs(hgt - INVALID_DATA_VALUE) < 1e-6*INVALID_DATA_VALUE:
            raise InvalidDataPointError('Invalid data point value [HGT]')
    elif reanalysis == "MERRA2" or reanalysis == "GEOS5":
        hgt = data[MERRA_HGT_PARMS[time]][layer]\
                [point.reanalysis_row-1][point.reanalysis_col-1]
        # For MERRA2 and GEOS5, invalid data is filled.
        if abs(hgt - INVALID_DATA_VALUE) < 1e-6*INVALID_DATA_VALUE:
            return INVALID_DATA_VALUE
    else:
        raise InvalidDataPointError('Unable to locate height value for {0}'.
                format(reanalysis))

    rad_lat = math.radians(point.lat)
    sin_lat = math.sin(rad_lat)
    cos_lat = math.cos(rad_lat)
    cos_2lat = math.cos(2.0 * rad_lat)

    radius = (1000.0 * math.sqrt(1.0 /
                                 (((cos_lat * cos_lat) *
                                   INV_R_MAX_SQRD) +
                                  ((sin_lat * sin_lat) *
                                   INV_R_MIN_SQRD))))

    gravity_ratio = ((9.80616 *
                      (1.0 -
                       (0.002637 * cos_2lat) +
                       (0.0000059 * (cos_2lat * cos_2lat)))) *
                     INV_STD_GRAVITY)

    return (hgt * radius) / (1000.0 * (gravity_ratio * radius - hgt))


MH_20 = 18.01534
MD_RY = 28.9644


def determine_rh_temp(data, point, layer, time, reanalysis):
    """Determines relative humidity and temperature for a layer in the point
       at time x

    Note: The layer value is the pressure, so needs conversion to float
          from int

    Args:
        data <dict>: Data structure for the parameters and pressure layers
        point <GridPointInfo>: The current point information
        layer <str>: The layer to get the relative humidity for
        time <int>: Time 0 or Time 1
        reanalysis <str>: Reanalysis source: NARR, MERRA2 or GEOS5

    Returns:
        <float>: Relative humidity for the point layer
        <float>: Temperature for the point layer
    """

    if reanalysis == "NARR":

        # Specific humidity at time
        spfh = data[NARR_SPFH_PARMS[time]][layer]\
                [point.reanalysis_row-1][point.reanalysis_col-1]
        # Temperature at time
        temp = data[NARR_TMP_PARMS[time]][layer]\
                [point.reanalysis_row-1][point.reanalysis_col-1]

        # For NARR, invalid data is an error case.
        if abs(spfh - INVALID_DATA_VALUE) < 1e-6*INVALID_DATA_VALUE:
            raise InvalidDataPointError('Invalid data point value [SPFH]')
        if abs(temp - INVALID_DATA_VALUE) < 1e-6*INVALID_DATA_VALUE:
            raise InvalidDataPointError('Invalid data point value [TMP]')

    elif reanalysis == "MERRA2" or reanalysis == "GEOS5":

        # Specific humidity at time
        spfh = data[MERRA_SPFH_PARMS[time]][layer]\
                [point.reanalysis_row-1][point.reanalysis_col-1]
        # Temperature at time
        temp = data[MERRA_TMP_PARMS[time]][layer]\
                [point.reanalysis_row-1][point.reanalysis_col-1]

        # For MERRA2 and GEOS5, invalid data is filled.
        if abs(spfh - INVALID_DATA_VALUE) < 1e-6*INVALID_DATA_VALUE:
            return (INVALID_DATA_VALUE, INVALID_DATA_VALUE)
        if abs(temp - INVALID_DATA_VALUE) < 1e-6*INVALID_DATA_VALUE:
            return (INVALID_DATA_VALUE, INVALID_DATA_VALUE)

    else:
        raise InvalidDataPointError('Unable to locate humidity value for {0}'.
                format(reanalysis))

    # Calculate vapor pressure at given temperature - hpa
    goff_pow_1 = math.pow(10.0, (11.344 * (1.0 - (temp / 373.16)))) - 1.0
    goff_pow_2 = math.pow(10.0, (-3.49149 * (373.16 / temp - 1.0))) - 1.0
    goff = (-7.90298 * (373.16 / temp - 1.0) +
            5.02808 * math.log10(373.16 / temp) -
            1.3816e-7 * goff_pow_1 +
            8.1328e-3 * goff_pow_2 +
            math.log10(1013.246))  # hPa

    # Calculate partial pressure
    ph_20 = ((spfh * float(layer) * MD_RY) /
             (MH_20 - spfh * MH_20 + spfh * MD_RY))

    # Calculate relative humidity
    relative_humidity = (ph_20 / pow(10.0, goff)) * 100.0

    return (relative_humidity, temp)


LayerInfo = namedtuple('LayerInfo', ('hgt', 'rh', 'temp'))


def interpolate_to_pressure_layers(data, point, layers, interp_factor,
                                   reanalysis):
    """Interpolate to each pressure layer for the point

    Args:
        data <dict>: Data structure for the parameters and pressure layers
        point <GridPointInfo>: The current point information
        layers [<str>]: All the pressure layers to load for each parameter
        interp_factor <float>: The interpolation factor to use
        reanalysis <str>: Reanalysis source: NARR, MERRA2 or GEOS5

    Returns:
        <dict>: dict[layers]->LayerInfo
    """

    interpolated_values = dict()
    for layer in layers:
        # Geometric height at t0
        hgt_m0 = determine_geometric_hgt(data, point, layer, 0, reanalysis)

        # Geometric height at t1
        hgt_m1 = determine_geometric_hgt(data, point, layer, 1, reanalysis)

        # Relative humidity at t0
        (rh_v0, temp_v0) = determine_rh_temp(data, point, layer, 0, reanalysis)

        # Relative humidity at t1
        (rh_v1, temp_v1) = determine_rh_temp(data, point, layer, 1, reanalysis)

        '''
        Linearly interpolate geometric height, relative humidity, and
        temperature for reanalysis points.  This is the reanalysis data
        corresponding to the acquisition date and scene center time of
        the Landsat data converted to appropriate values for MODTRAN
        runs.
        '''

        # Geometric height at acquisition date and scene center time
        if hgt_m0 == INVALID_DATA_VALUE or hgt_m1 == INVALID_DATA_VALUE:
            hgt = INVALID_DATA_VALUE
        else:
            hgt = (hgt_m0 + ((hgt_m1 - hgt_m0) * interp_factor))

        # Relative humidity at acquisition date and scene center time
        if rh_v0 == INVALID_DATA_VALUE or rh_v1 == INVALID_DATA_VALUE:
            relative_humidity = INVALID_DATA_VALUE
        else:
            relative_humidity = (rh_v0 + ((rh_v1 - rh_v0) * interp_factor))

        # Temperature at acquisition date and scene center time
        if temp_v0 == INVALID_DATA_VALUE or temp_v1 == INVALID_DATA_VALUE:
            temp = INVALID_DATA_VALUE
        else:
            temp = (temp_v0 + ((temp_v1 - temp_v0) * interp_factor))

        interpolated_values[layer] = LayerInfo(hgt=hgt, rh=relative_humidity,
                                               temp=temp)

    return interpolated_values


def get_latitude_longitude_strings(point):
    """Determine string versions of the latitude and longitude

    MODTRAN uses longitudinal degree values from 0 to 360 starting at
    Greenwich and moving west.  The following logic here fixes the
    longitude to be for MODTRAN.

    Args:
        point <GridPointInfo>: The current point information

    Returns:
        latitude <str>: The latitude in string format
        longitude <str>: The longitude in string format
    """

    latitude = '{0:06.3f}'.format(point.lat)
    longitude = None
    if point.lon < 0:
        '''
        We are a west longitude value so negate it to the positive
        equivalent value.
        i.e.   W40 normally represented as -40 would be changed to a
               positive 40 value.
        '''
        longitude = '{0:06.3f}'.format(-point.lon)
    else:
        '''
        We are a east longitude value so fix it to be greater than
        180 west.
        i.e.   E10 would be turned into W350, and be positive 350 not
               negative.
        '''
        longitude = '{0:06.3f}'.format(360.0 - point.lon)

    return (latitude, longitude)


def determine_layers_below_above(layers, values, elevation):
    """Determine the layers below and above the current elevation

    Args:
        layers [<str>]: All the pressure layers to load for each parameter
        values [<LayerInfo>]: All the interpolated layers information
        elevation <float>: The elevation to determine the layers for

    Returns
        <str>: Below layer ID
        <str>: Above layer ID
        <int>: Index to the above layer
    """

    index = 0
    index_below = 0
    index_above = len(layers) - 1

    for layer in layers:
        if values[layer].hgt >= elevation:
            index_below = index - 1
            index_above = index
            break

        index += 1

    # Only need to check for the low height condition
    # We should never have a height above our highest elevation
    if index_below < 0:
        index_below = 0
        index_above = 1

    return (layers[index_below], layers[index_above], index_above)


def interpolate_to_elevation(values, elevation, below, above):
    """Interpolate to our current elevation between the below and above layers

    Args:
        values [<LayerInfo>]: All the interpolated layers information
        elevation <float>: The current elevation
        below <str>: Below layer ID
        above <str>: Above layer ID

    Returns:
        <float>: Pressure for the elevation
        <float>: Temperature for the elevation
        <float>: Relative Humidity for the elevation
    """

    # Setup for linear interpolation between the heights
    hgt_interp_factor = \
        determine_interp_factor(elevation,
                                values[below].hgt,
                                values[above].hgt)

    # Linear interpolate pressure, temperature, and relative
    # humidity to elevation for lowest layer
    pressure = (float(below) +
                ((float(above) - float(below)) * hgt_interp_factor))

    temp = (values[below].temp +
            ((values[above].temp - values[below].temp) * hgt_interp_factor))

    rel_hum = (values[below].rh +
               ((values[above].rh - values[below].rh) * hgt_interp_factor))

    return (pressure, temp, rel_hum)


PressureLayerInfo = namedtuple('PressureLayerInfo',
                               ('hgt', 'pressure', 'temp', 'rh'))


def determine_base_layers_for_elev(layers, values, elevation):
    """Determine all of the base layers to use for the elevation

    Args:
        layers [<str>]: All the pressure layers
        values [<LayerInfo>]: All the interpolated layers information
        elevation <float>: The current elevation

    Returns:
        [PressureLayerInfo]: List of pressure layer information
    """

    # Determine layers above and below current elevation
    (below,
     above,
     index_above) = determine_layers_below_above(layers=layers,
                                                 values=values,
                                                 elevation=elevation)

    b_layers = list()

    '''
    MODTRAN throws an error when there are two identical layers in
    the tape5 file, if the current ground altitude and the next
    highest layer are close enough, don't add interpolated layer
    '''
    if math.fabs(elevation - values[above].hgt) > 0.001:
        # Not too close so generate and add it as the first element
        (pressure,
         temp,
         relative_humidity) = interpolate_to_elevation(values=values,
                                                       elevation=elevation,
                                                       below=below,
                                                       above=above)

        b_layers.append(PressureLayerInfo(hgt=elevation,
                                          pressure=pressure,
                                          temp=temp,
                                          rh=relative_humidity))

    # Remaining elements are the pressure layers above so add them
    for layer in layers[index_above:]:

        b_layers.append(PressureLayerInfo(hgt=values[layer].hgt,
                                          pressure=float(layer),
                                          temp=values[layer].temp,
                                          rh=values[layer].rh))

    return b_layers


def determine_all_layers_for_elev(std_atmos, layers, values, elevation):
    """Determine all of the layers to use for the elevation

    Args:
        std_atmos [StdAtmosInfo]: The standard atmosphere
        layers [<str>]: All the pressure layers
        values [<LayerInfo>]: All the interpolated layers information
        elevation <float>: The current elevation

    Returns:
        [PressureLayerInfo]: List of pressure layer information
    """

    s_layers = determine_base_layers_for_elev(layers=layers,
                                              values=values,
                                              elevation=elevation)

    # Determine maximum height of base layers and where the
    # standard atmosphere is greater than this
    first_index = None
    layer_index = 0
    for layer in std_atmos:
        # Check against the top base pressure layer
        if layer.hgt > s_layers[-1].hgt:
            first_index = layer_index
            break

        layer_index += 1

    second_index = first_index + 2

    '''
    If there are more than 2 layers above the highest reanalysis layer,
    then we need to interpolate a value between the highest reanalysis
    layer and the 2nd standard atmosphere layer above the reanalysis
    layers to create a smooth transition between the reanalysis layers
    and the standard upper atmosphere
    '''

    # Add an interpolated std layer
    if len(std_atmos[first_index:]) >= 3:

        # Setup for linear interpolation between the layers
        hgt = ((std_atmos[second_index].hgt + s_layers[-1].hgt) / 2.0)
        std_interp_factor = \
            determine_interp_factor(hgt, s_layers[-1].hgt,
                                    std_atmos[second_index].hgt)

        pressure = (s_layers[-1].pressure +
                    ((std_atmos[second_index].pressure -
                      s_layers[-1].pressure) * std_interp_factor))

        temp = (s_layers[-1].temp +
                ((std_atmos[second_index].temp -
                  s_layers[-1].temp) * std_interp_factor))

        rel_hum = (s_layers[-1].rh +
                   ((std_atmos[second_index].rh -
                     s_layers[-1].rh) * std_interp_factor))

        s_layers.append(PressureLayerInfo(hgt=hgt,
                                          pressure=pressure,
                                          temp=temp,
                                          rh=rel_hum))

    # Add the remaining standard atmosphere layers
    for layer in std_atmos[second_index:]:

        s_layers.append(PressureLayerInfo(hgt=layer.hgt,
                                          pressure=layer.pressure,
                                          temp=layer.temp,
                                          rh=layer.rh))

    return s_layers


def write_tape5_file(filename, head_data, body_data, tail_data):
    """Write the tape5 file with the specified information

    Args:
        filename <str>: The path and filename to create
        head_data <str>: The head information for the file
        body_data <str>: The body information for the file
        tail_data <str>: The tail information for the file
    """

    logger = logging.getLogger(__name__)

    logger.debug('Creating [{}]'.format(filename))

    with open(filename, 'w') as tape5_fd:
        tape5_fd.write(head_data)
        tape5_fd.write(body_data)
        tape5_fd.write(tail_data)


def generate_for_temp_alb_pairs(hgt_path, temp_head_data, body_data,
                                tail_data):
    """Generates all of the tape5 files for the temperature and albedo pairs

    Args:
        hgt_path <str>: Path to the height directory to place the temperature
                        and albedo pairs
        temp_head_data <str>: Template head data for the file
        body_data <str>: The body information for the file
        tail_data <str>: The tail information for the file
    """

    # Iterate through all the (temperature,albedo) pairs at which to
    # run MODTRAN and create the final required tape5 file for
    # each one
    for (temperature, albedo) in TEMP_ALBEDO_PAIRS:
        # Append the temperature and albedo to the path
        ta_path = os.path.join(hgt_path, temperature, albedo)

        # Now create before writing the tape5 file
        util.System.create_directory(ta_path)

        # Update the head section for the MODTRAN tape5 file with
        # current information
        head_data = temp_head_data.replace('tmp', temperature)
        head_data = head_data.replace('alb', albedo)

        write_tape5_file(filename=os.path.join(ta_path, TAPE5),
                         head_data=head_data,
                         body_data=body_data,
                         tail_data=tail_data)


def generate_tape5_for_elevation(std_atmos, layers, head_template, tail_data,
                                 point_path, values, elevation):
    """Generates all of the tape5 files for the elevation

    Args:
        std_atmos [StdAtmosInfo]: The standard atmosphere
        layers [<str>]: All the pressure layers
        head_template <str>: Template for the head of a tape5 file
        tail_data <str>: Updated tail data for a tape5 file
        point_path <str>: Path to the point directory to contain the tape5
                          files
        values [<LayerInfo>]: All the interpolated layers information
        elevation <float>: The current elevation
    """

    # Append the height directory
    hgt_path = os.path.join(point_path, '{0:05.3f}'.format(elevation))

    s_layers = determine_all_layers_for_elev(std_atmos=std_atmos,
                                             layers=layers,
                                             values=values,
                                             elevation=elevation)

    # Update the middle section for the MODTRAN tape5 file with
    # current information
    body_data = ''.join(['{0:10.3f}{1:10.3e}{2:10.3e}{3:10.3e}{4:10.3e}'
                         '{5:10.3e}{6:16s}\n'.format(x.hgt,
                                                     x.pressure,
                                                     x.temp,
                                                     x.rh,
                                                     0.0, 0.0,
                                                     'AAH             ')
                         for x in s_layers])

    # Update the head section for the MODTRAN tape5 file with
    # current information
    temp_head_data = head_template.replace('nml', str(len(s_layers)))
    temp_head_data = temp_head_data.replace('gdalt',
                                            '{0:05.3f}'
                                            .format(elevation))

    generate_for_temp_alb_pairs(hgt_path=hgt_path,
                                temp_head_data=temp_head_data,
                                body_data=body_data,
                                tail_data=tail_data)


def determine_elevations(grid_elevation_file, elevations, height):
    """Determine a list of adjusted elevations to use

    Args:
        grid_elevation_file <file>: File where 0 elevations are to be written
        elevations [<float>]: Initial list of elevations
        height <float>: Height for the bottom pressure layer

    Returns
        [<float>]: List of elevations to use
    """

    # Copy the elevations to iterate over
    new_elevations = [x for x in elevations]

    # Adjust the first element to be the geometric height unless it
    # is negative, then use 0.0
    if height < 0.0:
        new_elevations[0] = 0.0
    else:
        new_elevations[0] = height

    # Write the first elevation to the elevation file.  Write it with more
    # precision for science calculations and less to exactly match the
    # directory name
    grid_elevation_file.write('{0:05.8f} {0:05.3f}\n'.format(new_elevations[0]))

    return new_elevations


def build_ground_altitudes(ground_altitudes, espa_metadata):

    """Build list of ground altitudes based on the standard ground altitudes
       and the scene.  Also, write these altitudes to a file for later use

    Args:
        ground_altitudes <file>: List of altitudes to process
        espa_metadata <espa.metadata>: The metadata information for the input
    """

    # Determine the minimum and maximum scene elevations, excluding fill
    data_type = np.int16
    for band in espa_metadata.xml_object.bands.band:
        if (band.get('product') == 'elevation' and
                band.get('category') == 'image'):
            filename = str(band.file_name)
            elevation_fill = int(band.get('fill_value'))

    input_data = np.fromfile(filename, dtype=data_type)
    masked_input_data = np.ma.masked_array(input_data, input_data == \
        elevation_fill)
    elevation_minimum = np.min(masked_input_data)
    elevation_maximum = np.max(masked_input_data)

    # If minimum or maximum is less than 0, make it 0
    if elevation_minimum < 0.0:
        elevation_minimum = 0.0
    if elevation_maximum < 0.0:
        elevation_maximum = 0.0

    # Scale the elevation results to match the reanalysis data (scale m to km)
    scaled_elevation_minimum = elevation_minimum * 0.001
    scaled_elevation_maximum = elevation_maximum * 0.001

    # Find the standard altitudes that bound the minimum and maximum elevations
    min_index = bisect(GROUND_ALT, scaled_elevation_minimum) - 1
    max_index = bisect(GROUND_ALT, scaled_elevation_maximum)

    # Make sure elevation indexes are within range (the actual elevations could
    # fall outside the standard range).  Also, we don't need to put the 0.0
    # elevation in the list twice, so avoid that for low elevation scenes
    num_elevations = len(GROUND_ALT)
    if min_index <= 0:
        min_index = 1
    if min_index >= num_elevations:
        min_index = num_elevations - 1
    if max_index <= 0:
        max_index = 1
    if max_index >= num_elevations:
        max_index = num_elevations - 1

    # Build the elevations list like the standard elevations list, but only
    # with elevations we need to process based on this scene's elevations.
    # We always need the 0 elevation
    ground_altitudes.append(0.0)
    for index in range(min_index, max_index + 1):
        ground_altitudes.append(GROUND_ALT[index])

    # Make a file that shows what elevations are actually being used
    with open(MODTRAN_ELEVATION_NAME, 'w') as modtran_elevation_file:
        modtran_elevation_file.write('{0}\n'.format(len(ground_altitudes)))
        for altitude in ground_altitudes:
            modtran_elevation_file.write('{0}\n'.format(altitude))
    modtran_elevation_file.close()


def least_squares_estimate(x, y, p):

    """Calculate new interpolated or extrapolated value for one of the 
       atmospheric parameters where there is missing data.

    Args:
        x <2d array>: Based on nearby pressure level values
        y <2d array>: Based on nearby parameter values at these pressure levels
        p <int>: Pressure layer of the value to be interpolated or extrapolated
    """

    mult1 = np.dot(x, x.conj().transpose())
    mult2 = np.dot(y, x.conj().transpose())
    inv_mult1 = np.linalg.inv(mult1)
    a = np.dot(mult2, inv_mult1)
    m = a[0,1]
    b = a[0,0]

    return (m * p + b)


def extrapolate_bottom(layer_index, values, field):

    """The bottom element of the selected atmospheric parameter is missing, so
       extrapolate its value from higher elements.
       Used only for MERRA2 or GEOS5.

    Args:
        layer_index <int>: Pressure layer index 
        values <dict>: Atmospheric parameters, dict[layers]->LayerInfo
        field <string>: The named tuple field for the parameter we are updating
    """

    # Extrapolate using the first 3 present values.  This assumes that once
    # you find the first present value, the next 2 are also present.
    level = layer_index + 1
    while getattr(values[MERRA_PRESSURE_LAYERS[level]], field) \
        == INVALID_DATA_VALUE: 
        level = level + 1

    x = np.ones((2,3))
    x[1] = MERRA_PRESSURE_LAYERS[level:level+3]

    y = np.zeros((1,3))
    y[0,0] = getattr(values[MERRA_PRESSURE_LAYERS[level]], field)
    y[0,1] = getattr(values[MERRA_PRESSURE_LAYERS[level+1]], field)
    y[0,2] = getattr(values[MERRA_PRESSURE_LAYERS[level+2]], field)

    # Pull out the pressure layer for readability.
    p = MERRA_PRESSURE_LAYERS[layer_index]

    new_value = least_squares_estimate(x, y, p)

    # Replace the named tuple with a new one with 1 value updated.
    values[p] = values[p]._replace(**{field : new_value})

    return values 


def extrapolate_top(layer_index, values, field):

    """One of the top elements of the selected atmospheric parameter is missing,
       so extrapolate its value from lower elements.
       Used only for MERRA2 or GEOS5.

    Args:
        layer_index <int>: Pressure layer index 
        values <dict>: Atmospheric parameters, dict[layers]->LayerInfo
        field <string>: The named tuple field for the parameter we are updating
    """

    # Extrapolate using 2 present values below.  Fill in values from bottom to 
    # top.  Values below should always be present.
    x = np.ones((2,2))
    x[1,0] = MERRA_PRESSURE_LAYERS[layer_index - 2]
    x[1,1] = MERRA_PRESSURE_LAYERS[layer_index - 1]

    y = np.zeros((1,2))
    y[0,0] = getattr(values[MERRA_PRESSURE_LAYERS[layer_index - 2]], field)
    y[0,1] = getattr(values[MERRA_PRESSURE_LAYERS[layer_index - 1]], field)

    # Pull out the pressure layer for readability.
    p = MERRA_PRESSURE_LAYERS[layer_index]

    new_value = least_squares_estimate(x, y, p)

    # Replace named tuple with new one with 1 value updated.
    values[p] = values[p]._replace(**{field : new_value})

    return values 


def interpolate_middle(above, layer_index, values, field):

    """A middle element of the selected atmospheric parameter is missing, so
       interpolate its value from lower and higher elements.
       Used only for MERRA2 or GEOS5.

    Args:
        above <int>: Index of closest element above this one with valid data 
        layer_index <int>: Pressure layer index 
        values <dict>: Atmospheric parameters, dict[layers]->LayerInfo
        field <string>: The named tuple field for the parameter we are updating
    """

    below = layer_index - 1

    x = np.ones((2,2))
    x[1,0] = MERRA_PRESSURE_LAYERS[below]
    x[1,1] = MERRA_PRESSURE_LAYERS[above]

    y = np.zeros((1,2))
    y[0,0] = getattr(values[MERRA_PRESSURE_LAYERS[below]], field)
    y[0,1] = getattr(values[MERRA_PRESSURE_LAYERS[above]], field)

    # Pull out the pressure layer for readability.
    p = MERRA_PRESSURE_LAYERS[layer_index]

    new_value = least_squares_estimate(x, y, p)

    # Replace named tuple with a new one with 1 value updated.
    values[p] = values[p]._replace(**{field : new_value})

    return values 


def fill_missing_values(layer_index, values, field):

    """For the selected atmospheric parameter, extrapolate or interpolate to 
       fill it a missing value.  In some cases fill in nearby values, too.
       Used only for MERRA2 or GEOS5.

    Args:
        layer_index <int>: Pressure layer index 
        values <dict>: Atmospheric parameters, dict[layers]->LayerInfo
        field <string>: The named tuple field for the parameter we are checking 
    """

    if layer_index == 0:

        # The first value is missing.
        values = extrapolate_bottom(layer_index, values, field)
    else:
        above = layer_index

        # Find the first present value above the missing value.
        while getattr(values[MERRA_PRESSURE_LAYERS[above]], field) \
            == INVALID_DATA_VALUE: 
            above = above + 1

            if above == len(MERRA_PRESSURE_LAYERS):

                # There are no present values above.
                values = extrapolate_top(layer_index, values, field)
                break

        # If there are present values above, interpolate. 
        if above != len(MERRA_PRESSURE_LAYERS):
            values = interpolate_middle(above, layer_index, values, field)

    return values


def generate_tape5_files_for_point(grid_elevation_file, std_atmos, data, point,
                                   interp_factor, doy_str, head_template,
                                   tail_template, ground_altitudes, reanalysis):
    """Generate tape5 file for the current point

    Args:
        grid_elevation_file <file>: File where 0 elevations are to be written
        std_atmos [StdAtmosInfo]: The standard atmosphere
        data <dict>: Data structure for the parameters and pressure layers
        point <GridPointInfo>: The current point information
        interp_factor <float>: The interpolation factor to use
        head_template <str>: The template for the head of the tape5 file
        tail_template <str>: The template for the tail of the tape5 file
        ground_altitudes [<float>]: The standard altitudes we need to process
        reanalysis <str>: Reanalysis source: NARR, MERRA2 or GEOS5
    """

    logger = logging.getLogger(__name__)

    # Define the point directory
    point_path = ('{0:03}_{1:03}_{2:03}_{3:03}'
                  .format(point.row, point.col,
                          point.reanalysis_row, point.reanalysis_col))

    (latitude, longitude) = get_latitude_longitude_strings(point=point)
    logger.debug('MODTRAN latitude [{}]'.format(latitude))
    logger.debug('MODTRAN longitude [{}]'.format(longitude))

    # Update the tail section for the MODTRAN tape5 file with current
    # information
    tail_data = tail_template.replace('latitu', latitude)
    tail_data = tail_data.replace('longit', longitude)
    tail_data = tail_data.replace('jay', doy_str)

    if reanalysis == "NARR":
        PRESSURE_LAYERS = NARR_PRESSURE_LAYERS
    elif reanalysis == "MERRA2" or reanalysis == "GEOS5":
        PRESSURE_LAYERS = MERRA_PRESSURE_LAYERS
    else:
        PRESSURE_LAYERS = []

    values = interpolate_to_pressure_layers(data=data,
                                            point=point,
                                            layers=PRESSURE_LAYERS,
                                            interp_factor=interp_factor,
                                            reanalysis=reanalysis)

    if reanalysis == "MERRA2" or reanalysis == "GEOS5":
        # Check for missing values at each pressure level and fill them.
        # This is only done for MERRA2 and GEOS5.  For NARR, missing values are
        # considered an error case.
        for (layer_index, layer) in enumerate(PRESSURE_LAYERS):

            # If a value is missing, extrapolate or interpolate to fill it. 
            if values[layer].temp == INVALID_DATA_VALUE:
                values = fill_missing_values(layer_index, values, "temp")
            if values[layer].hgt == INVALID_DATA_VALUE:
                values = fill_missing_values(layer_index, values, "hgt")
            if values[layer].rh == INVALID_DATA_VALUE:
                values = fill_missing_values(layer_index, values, "rh")

    elevations = determine_elevations(grid_elevation_file,
                                      elevations=ground_altitudes,
                                      height=values[PRESSURE_LAYERS[0]].hgt)

    for elevation in elevations:
        generate_tape5_for_elevation(std_atmos=std_atmos,
                                     layers=PRESSURE_LAYERS,
                                     head_template=head_template,
                                     tail_data=tail_data,
                                     point_path=point_path,
                                     values=values,
                                     elevation=elevation)


def generate_modtran_tape5_files(espa_metadata, data_path, std_atmos,
                                 grid_points, reanalysis):
    """
    Args:
        espa_metadata <espa.metadata>: The metadata information for the input
        data_path <str>: The directory for the reanalysis data files
        std_atmos [StdAtmosInfo]: The standard atmosphere
        grid_points [GridPointInfo]: List of the grid point information
        reanalysis <str>: Reanalysis source: NARR, MERRA2 or GEOS5
    """

    # Load the MODTRAN head and tail template files
    head_template = load_modtran_template_file(data_path, MODTRAN_HEAD)
    tail_template = load_modtran_template_file(data_path, MODTRAN_TAIL)

    # Load the pressure layers into a data structure
    if reanalysis == "NARR":
        data = load_pressure_layers(parameters=NARR_PARAMETERS,
                                    layers=NARR_PRESSURE_LAYERS,
                                    reanalysis=reanalysis)
    elif reanalysis == "MERRA2" or reanalysis == "GEOS5":
        data = load_pressure_layers(parameters=MERRA_PARAMETERS,
                                    layers=MERRA_PRESSURE_LAYERS,
                                    reanalysis=reanalysis)
    else:
        data = None

    # Get the dates and determine the day-of-year
    (acq_date, t0_date, t1_date) = util.REANALYSIS.dates(espa_metadata)
    doy_str = '{0:03d}'.format(acq_date.timetuple().tm_yday)

    # Setup for linear interpolation between the dates in seconds
    interp_factor = (float((acq_date - t0_date).seconds) /
                     float((t1_date - t0_date).seconds))

    # Build the list of altitudes to process based on scene elevations
    ground_alts = []
    build_ground_altitudes(ground_alts, espa_metadata)

    with open(GRID_ELEVATION_NAME, 'w') as grid_elevation_file:
        for point in grid_points:
            if point.run_modtran:
                generate_tape5_files_for_point(grid_elevation_file,
                                               std_atmos=std_atmos,
                                               data=data,
                                               point=point,
                                               interp_factor=interp_factor,
                                               doy_str=doy_str,
                                               head_template=head_template,
                                               tail_template=tail_template,
                                               ground_altitudes=ground_alts,
                                               reanalysis=reanalysis)
    grid_elevation_file.close()


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

    logger.info('*** Begin MODTRAN Tape5 Generation ***')

    # XML Metadata
    espa_metadata = Metadata()
    espa_metadata.parse(xml_filename=args.xml_filename)

    # Load the grid information
    (grid_points, dummy1, dummy2) = read_grid_points()

    # Load the standard atmospheric layers information
    std_atmos = [layer for layer in
                 load_std_atmosphere(data_path=args.data_path)]

    generate_modtran_tape5_files(espa_metadata=espa_metadata,
                                 data_path=args.data_path,
                                 std_atmos=std_atmos,
                                 grid_points=grid_points,
                                 reanalysis=args.reanalysis)

    logger.info('*** MODTRAN Tape5 Generation - Complete ***')


if __name__ == '__main__':
    main()
