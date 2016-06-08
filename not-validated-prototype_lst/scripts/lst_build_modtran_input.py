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
import math
import array
from argparse import ArgumentParser
from collections import namedtuple
from osgeo import gdal, osr

from espa import Metadata
from lst_exceptions import InvalidParameterFileError
import lst_utilities as util

from lst_grid_points import read_grid_points


NARR_ROWS = 277
NARR_COLS = 349

PRESSURE_LAYERS = [1000, 975, 950, 925, 900,
                   875, 850, 825, 800, 775,
                   750, 725, 700, 650, 600,
                   550, 500, 450, 400, 350,
                   300, 275, 250, 225, 200,
                   175, 150, 125, 100]

'''
TODO TODO TODO TODO TODO TODO TODO TODO
TODO TODO TODO TODO TODO TODO TODO TODO

I think if we know more information about the data points that will need a
NARR point, the number of ALTITUDES needed to run through MODTRAN for a point
could be reduced, further reducing the number of MODTRAN runs.

I think it would be significant for flat lands to only have to run 6 MODTRAN
runs for a point instead of 27.

Mountainous land could end up running all.

I think I have the information and can add it to the grid points data.

TODO TODO TODO TODO TODO TODO TODO TODO
TODO TODO TODO TODO TODO TODO TODO TODO
'''
GROUND_ALT = [0.0, 0.6, 1.1,
              1.6, 2.1, 2.6,
              3.1, 3.6, 4.05]

# We can use strings for these in this code to make things simpler for
# directory creation
TEMPERATURES = ['273', '310', '000']
ALBEDOS = ['0.0', '0.0', '0.1']
TEMP_ALBEDO_PAIRS = zip(TEMPERATURES, ALBEDOS)

HGT_PARMS = ['HGT_t0', 'HGT_t1']
SPFH_PARMS = ['SPFH_t0', 'SPFH_t1']
TMP_PARMS = ['TMP_t0', 'TMP_t1']
PARAMETERS = HGT_PARMS + SPFH_PARMS + TMP_PARMS

MODTRAN_HEAD = 'modtran_head.txt'
MODTRAN_MIDDLE = 'modtran_middle.txt'
MODTRAN_TAIL = 'modtran_tail.txt'


def retrieve_command_line_arguments():
    """Read arguments from the command line

    Returns:
        args <arguments>: The arguments read from the command line
    """

    parser = ArgumentParser(description='Builds MODTRAN input data files for'
                                        ' a pre-determined set of points')

    parser.add_argument('--xml',
                        action='store', dest='xml_filename',
                        required=False, default=None,
                        help='The XML metadata file to use')

    parser.add_argument('--lst_data_dir',
                        action='store', dest='lst_data_dir',
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
        print(util.Version.version_text())
        sys.exit(0)  # EXIT SUCCESS

    if args.xml_filename is None:
        raise Exception('--xml must be specified on the command line')

    if args.xml_filename == '':
        raise Exception('The XML metadata filename provided was empty')

    if args.lst_data_dir is None:
        raise Exception('--lst_data_dir must be specified on the command line')

    if args.lst_data_dir == '':
        raise Exception('The LST data directory provided was empty')

    return args


def read_narr_pressure_file(parameter, layer):
    """
    Args:
        TODO TODO TODO
    """

    logger = logging.getLogger(__name__)

    filename = os.path.join(parameter, '.'.join([str(layer), 'txt']))
    logger.debug('Reading NARR Pressure File [{}]'.format(filename))

    point_values = None
    with open(filename, 'r') as data_fd:
        point_values = [[float(data_fd.readline())
                         for row in xrange(NARR_ROWS)]
                        for col in xrange(NARR_COLS)]

    return point_values


STD_ATMOS_FILENAME = 'std_mid_lat_summer_atmos.txt'
StdAtmosInfo = namedtuple('StdAtmosInfo',
                          ('hgt', 'pressure', 'temp', 'rh'))


def read_std_mid_lat_summer_atmos_file(lst_data_dir):
    """
    Args:
        TODO TODO TODO
    """

    logger = logging.getLogger(__name__)

    filename = os.path.join(lst_data_dir, STD_ATMOS_FILENAME)
    logger.debug('Reading Standard Atmosphere File [{}]'.format(filename))

    with open(filename, 'r') as data_fd:
        for line in data_fd.readlines():
            (hgt, pressure, temp, rh) = line.strip().split()
            yield(StdAtmosInfo(hgt=float(hgt), pressure=float(pressure),
                               temp=float(temp), rh=float(rh)))


EQUATORIAL_RADIUS = 6378137.0
POLAR_RADIUS = 6356752.3142
STANRDARD_GRAVITY_IN_M_PER_SEC_SQRD = 9.80665
EQUATORIAL_RADIUS_IN_KM = EQUATORIAL_RADIUS / 1000.0
POLAR_RADIUS_IN_KM = POLAR_RADIUS / 1000.0
INV_R_MIN_SQRD = 1.0 / (POLAR_RADIUS_IN_KM * POLAR_RADIUS_IN_KM)
INV_R_MAX_SQRD = 1.0 / (EQUATORIAL_RADIUS_IN_KM * EQUATORIAL_RADIUS_IN_KM)
INV_STD_GRAVITY = 1.0 / STANRDARD_GRAVITY_IN_M_PER_SEC_SQRD


def determine_geometric_height(data, point, layer, time):
    """Converts from geopotential to geometric height

    Geopotential height is converted to geometrix height for both t0 and t1
    heights

    Args:
        TODO TODO TODO

    Returns:
        TODO TODO TODO
    """

    # Geopotential height at time
    hgt = data[HGT_PARMS[time]][layer][point.narr_col][point.narr_row]

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


def determine_relative_humidities(data, point, layer, time):
    """

    Note: The layer value is the pressure, so needs conversion to float

    Args:
        TODO TODO TODO

    Returns:
        TODO TODO TODO
    """

    # Specific humidity at time
    spfh = data[SPFH_PARMS[time]][layer][point.narr_col][point.narr_row]
    # Temperature at time
    tmp = data[TMP_PARMS[time]][layer][point.narr_col][point.narr_row]

    # calculate vapor pressure at given temperature - hpa
    goff_pow_1 = math.pow(10.0, (11.344 * (1.0 - (tmp / 373.16)))) - 1.0
    goff_pow_2 = math.pow(10.0, (-3.49149 * (373.16 / tmp - 1.0))) - 1.0
    goff = (-7.90298 * (373.16 / tmp - 1.0) +
            5.02808 * math.log10(373.16 / tmp) -
            1.3816e-7 * goff_pow_1 +
            8.1328e-3 * goff_pow_2 +
            math.log10(1013.246))  # hPa

    # calculate partial pressure
    ph_20 = ((spfh * float(layer) * MD_RY) /
             (MH_20 - spfh * MH_20 + spfh * MD_RY))

    # calculate relative humidity
    return (ph_20 / pow(10.0, goff)) * 100.0


def get_latitude_longitude_strings(point):
    """Determine string versions of the latitude and longitude

    MODTRAN uses longitudinal degree values from 0 to 360 starting at
    Greenwich and moving west.  The following logic here fixes the
    longitude to be for MODTRAN.

    Args:
        TODO TODO TODO

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


def butter:
    pass


def generate_modtran_tape5_files(espa_metadata, lst_data_dir, std_atmosphere,
                                 grid_points, grid_rows, grid_cols):
    """
    Args:
        TODO TODO TODO
    """

    logger = logging.getLogger(__name__)

    # Load the head and tail template files
    filename = os.path.join(lst_data_dir, MODTRAN_HEAD)
    head_template = None
    with open(filename, 'r') as template_fd:
        head_template = template_fd.read()

    filename = os.path.join(lst_data_dir, MODTRAN_TAIL)
    tail_template = None
    with open(filename, 'r') as template_fd:
        tail_template = template_fd.read()

    # Read the NARR pressure layers into a data structure
    data = dict()
    for parameter in PARAMETERS:
        data[parameter] = dict()
        for layer in PRESSURE_LAYERS:
            data[parameter][layer] = dict()
            data[parameter][layer] = read_narr_pressure_file(parameter, layer)

    (acq_date, t0_date, t1_date) = util.NARR.dates(espa_metadata)
    doy_str = '{0:03d}'.format(acq_date.timetuple().tm_yday)

    # Setup for linear interpolation between the dates in seconds
    # Determine divisor in seconds
    inv_t1_minus_t0 = 1.0 / (t1_date - t0_date).seconds
    # Determine difference to min time
    ta_minus_t0 = (acq_date - t0_date).seconds

    interp_factor = ta_minus_t0 * inv_t1_minus_t0

    # Determine geometric height and relative humidity
    for point in grid_points:
        if point.run_modtran:
            # Define the point directory
            point_path = ('{0:03}_{1:03}_{2:03}_{3:03}'
                          .format(point.row, point.col,
                                  point.narr_row, point.narr_col))

            (latitude, longitude) = get_latitude_longitude_strings(point)
            logger.debug('MODTRAN latitude [{}]'.format(latitude))
            logger.debug('MODTRAN longitude [{}]'.format(longitude))

            # Update the tail section for the MODTRAN tape5 file with current
            # information
            tail_data = tail_template.replace('latitu', latitude)
            tail_data = tail_data.replace('longit', longitude)
            tail_data = tail_data.replace('jay', doy_str)

            hgt_v = dict()
            rh_v = dict()
            tmp_v = dict()
            for layer in PRESSURE_LAYERS:
                # Geometric height at t0
                hgt_m0 = determine_geometric_height(data, point, layer, 0)
                # print(hgt_m0)

                # Geometric height at t1
                hgt_m1 = determine_geometric_height(data, point, layer, 1)
                # print(hgt_m1)

                # Relative humitidy at t0
                rh_v0 = determine_relative_humidities(data, point, layer, 0)
                # print(rh_v0)

                # Relative humitidy at t1
                rh_v1 = determine_relative_humidities(data, point, layer, 0)
                # print(rh_v1)

                '''
                Linearly interpolate geometric height, relative humidity, and
                temperature for NARR points.  This is the NARR data
                corresponding to the acquisition date and scene center time of
                the Landsat data converted to appropriate values for MODTRAN
                runs.
                '''

                # Geometric height at acquisition date and scene center time
                hgt_v[layer] = (hgt_m0 + ((hgt_m1 - hgt_m0) * interp_factor))

                # Relative humidity at acquisition date and scene center time
                rh_v[layer] = (rh_v0 + ((rh_v1 - rh_v0) * interp_factor))

                # Temperature at acquisition date and scene center time
                tmp_v[layer] = (tmp_v0 + ((tmp_v1 - tmp_v0) * interp_factor))

            # Copy the elevations to iterate over
            elevations = [x for x in GROUND_ALT]
            # Adjust the first element to be the geometric height unless it
            # is negative, then use 0.0
            if hgt_v[PRESSURE_LAYERS[0]] < 0:
                elevations[0] = 0.0
            else:
                elevations[0] = hgt_v[PRESSURE_LAYERS[0]]

            for elevation in elevations:
                # Append the height directory
                hgt_path = os.path.join(point_path,
                                        '{0:05.3f}'.format(elevation))

                # Determine layers below current elevation and closest
                # layer above and below
                index = 0
                for layer in PRESSURE_LAYERS:
                    if hgt_v[layer] >= elevation:
                        index_below = index - 1
                        index_above = index

                    index += 1

                # Only need to check for the low height condition
                # We should never have a height above our highest elevation
                if index_below < 0:
                    index_below = 0
                    index_above = 1

                layer_below = PRESSURE_LAYERS[index_below]
                layer_above = PRESSURE_LAYERS[index_above]

                # Setup for linear interpolation between the heights
                inv_ha_minus_hb = 1.0 / (hgt_v[layer_above] -
                                         hgt_v[layer_below])
                e_minus_hb = elevation - hgt_v[layer_below]
                hgt_interp_factor = e_minus_hb * inv_ha_minus_hb

                # Linear interpolate pressure, temperature, and relative
                # humidity to elevation for lowest layer
                new_pressure = (float(layer_below) +
                                ((float(layer_above) - float(layer_below)) *
                                 hgt_interp_factor))
                new_temp = (tmp_v[layer_below] +
                            ((tmp_v[layer_above] - tmp_v[layer_below]) *
                             hgt_interp_factor))
                new_rh = (rh_v[layer_below] +
                          ((rh_v[layer_above] - rh_v[layer_below]) *
                           hgt_interp_factor))

                # Create arrays of values containing only the layers to be
                # included in current tape5 file
                temp_hgt = list()
                temp_pressure = list()
                temp_temp = list()
                temp_rh = list()

                # First element is the recently figured-out values
                temp_hgt.append(elevation)
                temp_pressure.append(new_pressure)
                temp_temp.append(new_temp)
                temp_rh.append(new_rh)

                # Remaining elements are the pressure layers above
                for layer in PRESSURE_LAYERS[index_above:]:
                    temp_hgt.append(hgt_v[layer])
                    temp_pressure.append(float(layer))
                    temp_temp.append(tmp_v[layer])
                    temp_rh.append(rh_v[layer])

                '''
                MODTRAN throws an error when there are two identical layers in
                the tape5 file, if the current ground altitude and the next
                highest layer are close enough, eliminate interpolated layer
                '''
                if math.fabs(elevation - hgt_v[layer_above]) < 0.001:
                    # Simply remove the first element of the array because
                    # that is the element we just tested
                    del temp_height[0]
                    del temp_pressure[0]
                    del temp_temp[0]
                    del temp_rh[0]

                # Determine maximum height of NARR layers and where the
                # standard atmosphere is greater than this
                first_index = None
                layer_index = 0
                for layer in std_atmosphere:
                    # Check against the top pressure layer we have saved
                    if layer.hgt > temp_hgt[-1]:
                        first_index = layer_index
                        break

                    layer_index += 1
                second_index = first_index + 1

                '''
                If there are more than 2 layers above the highest NARR layer,
                then we need to interpolate a value between the highest NARR
                layer and the 2nd standard atmosphere layer above the NARR
                layers to create a smooth transition between the NARR layers
                and the standard upper atmosphere
                '''
                if len(std_atmosphere[first_index:]) >= 3:

                    # Setup for linear interpolation between the layers
                    inv_s_minus_l = 1.0 / (std_atmosphere[second_index].hgt -
                                           temp_hgt[-1])
                    std_height = ((std_atmosphere[second_index].hgt +
                                   temp_hgt[-1]) / 2.0)
                    h_minus_l = std_height - temp_hgt[-1]

                    std_interp_factor = h_minus_l * inv_s_minus_l

                    std_pressure = (temp_pressure[-1] +
                                    ((std_atmosphere[second_index].pressure -
                                      temp_pressure[-1]) * std_interp_factor))

                    std_temp = (temp_temp[-1] +
                                ((std_atmosphere[second_index].temp -
                                  temp_temp[-1]) * std_interp_factor))

                    std_rh = (temp_rh[-1] +
                              ((std_atmosphere[second_index].rh -
                                temp_rh[-1]) * std_interp_factor))

                    temp_hgt.append(std_height)
                    temp_pressure.append(std_pressure)
                    temp_temp.append(std_temp)
                    temp_rh.append(std_rh)

                # Add the remaining standard atmosphere layers
                for layer in std_atmosphere[second_index:]:
                    temp_hgt.append(layer.hgt)
                    temp_pressure.append(layer.pressure)
                    temp_temp.append(layer.temp)
                    temp_rh.append(layer.rh)

                # Update the middle section for the MODTRAN tape5 file with
                # current information
                num_layers = len(temp_hgt)
                middle_data = ''
                for index in xrange(num_layers):
                    middle_data = ''.join([middle_data,
                                           '{0:10.3f}'
                                           '{1:10.3e}'
                                           '{2:10.3e}'
                                           '{3:10.3e}'
                                           '{4:10.3e}'
                                           '{5:10.3e}'
                                           '{6:16s}\n'
                                           .format(temp_hgt[index],
                                                   temp_pressure[index],
                                                   temp_temp[index],
                                                   temp_rh[index],
                                                   0.0, 0.0,
                                                   'AAH             ')])

                # Update the head section for the MODTRAN tape5 file with
                # current information
                temp_head_data = head_template.replace('nml', str(num_layers))
                temp_head_data = temp_head_data.replace('gdalt',
                                                        '{0:05.3f}'
                                                        .format(elevation))

                # Iterate through all the temperature,albedo pairs at which to
                # run MODTRAN and create the final required tape5 file for
                # each one
                for (temperature, albedo) in TEMP_ALBEDO_PAIRS:
                    # Append the temperature and albedo to the path
                    ta_path = os.path.join(hgt_path, temperature, albedo)
                    # Now create before writing the tape5 file
                    logger.info('Creating [{}]'.format(ta_path))
                    util.System.create_directory(ta_path)

                    # Update the head section for the MODTRAN tape5 file with
                    # current information
                    head_data = temp_head_data.replace('nml', str(num_layers))
                    head_data = head_data.replace('tmp', temperature)
                    head_data = head_data.replace('alb', albedo)

                    with open('tape5', 'w') as tape5_fd:
                        tape5_fd.write(head_data)
                        tape5_fd.write(middle_data)
                        tape5_fd.write(tail_data)


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

    # XML Metadata
    espa_metadata = Metadata()
    espa_metadata.parse(xml_filename=args.xml_filename)

    # Load the grid information
    (grid_points, grid_rows, grid_cols) = read_grid_points()

    # Load the standard atmospheric layers information
    std_atmosphere = [layer for layer in
                      read_std_mid_lat_summer_atmos_file(args.lst_data_dir)]

    generate_modtran_tape5_files(espa_metadata, args.lst_data_dir,
                                 std_atmosphere=std_atmosphere,
                                 grid_points=grid_points,
                                 grid_rows=grid_rows,
                                 grid_cols=grid_cols)


if __name__ == '__main__':
    main()
