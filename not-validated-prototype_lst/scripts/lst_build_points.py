#! /usr/bin/env python

'''
    File: lst_build_points.py

    Purpose: Builds a directory structure of points and required information
             to be used for generation and interpolation of atmospheric
             information by follow on applications.

    Project: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    License: NASA Open Source Agreement 1.3
'''

import os
import sys
import logging
import math
import numpy as np
from argparse import ArgumentParser
from collections import namedtuple
from osgeo import gdal, osr

from espa import Metadata
from lst_exceptions import MissingBandError
import lst_utilities as util


NARR_COORDINATES_FILENAME = 'narr_coordinates.txt'
RADIUS_EARTH_IN_METERS = 6378137.0
DIAMETER_EARTH_IN_METERS = RADIUS_EARTH_IN_METERS * 2


GdalInfo = namedtuple('GdalInfo',
                      ('data_ds', 'data_srs',
                       'data_transform', 'll_srs',
                       'll_to_data', 'data_to_ll',
                       'nlines', 'nsamps', 'fill_value'))
DataBoundInfo = namedtuple('DataBoundInfo',
                           ('north_lat', 'south_lat', 'east_lon', 'west_lon'))
CoordinateInfo = namedtuple('CoordinateInfo', ('col', 'row'))
MinMaxRowColInfo = namedtuple('MinMaxRowColInfo',
                              ('min_row', 'max_row', 'min_col', 'max_col'))
PointInfo = namedtuple('PointInfo',
                       ('col', 'row',
                        'latitude', 'longitude',
                        'map_y', 'map_x'))


def retrieve_command_line_arguments():
    """Read arguments from the command line

    Returns:
        args <arguments>: The arguments read from the command line
    """

    # Create a command line arugement parser
    parser = ArgumentParser(description='Builds a directory structure of'
                                        ' points and required information'
                                        ' suitable for MODTRAN execution.')

    # ---- Add Arguments ----
    parser.add_argument('--xml',
                        action='store', dest='xml_filename',
                        required=False, default=None,
                        help='The XML metadata file to use')

    parser.add_argument('--lst_data_dir',
                        action='store', dest='lst_data_dir',
                        required=False, default=None,
                        help='Output debug messages and/or keep debug data')

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

    # Verify that the --lst_data_dir parameter was specified
    if args.lst_data_dir is None:
        raise Exception('--lst_data_dir must be specified on the'
                        ' command line')

    # Verify that the LST data directory provided is not an empty string
    if args.lst_data_dir == '':
        raise Exception('The LST data directory provided was empty')

    return args


def initialize_gdal_objects(espa_metadata):
    """Initialize the gdal objects and information we will need later
    """

    logger = logging.getLogger(__name__)

    # Use the brightness temperature band as the source for size and
    # projection information
    nlines = None
    nsamps = None
    fill_value = None
    toa_bt_filename = None
    for band in espa_metadata.xml_object.bands.band:
        if (band.get('product') == 'toa_bt' and
                band.get('category') == 'image'):

            toa_bt_filename = str(band.file_name)
            nlines = int(band.get('nlines'))
            nsamps = int(band.get('nsamps'))
            fill_value = int(band.get('fill_value'))
            break

    if toa_bt_filename is None:
        raise MissingBandError('Missing TOA Brightness Temperature Band')

    logger.debug('Filename: {}'.format(toa_bt_filename))

    data_ds = gdal.Open(toa_bt_filename)

    # Create the data SRS
    data_srs = osr.SpatialReference()
    data_srs.ImportFromWkt(data_ds.GetProjection())
    logger.debug('Data WKT: {}'.format(data_ds.GetProjection()))

    # Create the Geographic SRS
    ll_srs = data_srs.CloneGeogCS()

    # Create the coordinate transformations
    ll_to_data = osr.CoordinateTransformation(ll_srs, data_srs)
    data_to_ll = osr.CoordinateTransformation(data_srs, ll_srs)

    # Get the data transform
    data_transform = data_ds.GetGeoTransform()

    return GdalInfo(data_ds=data_ds, data_srs=data_srs,
                    data_transform=data_transform, ll_srs=ll_srs,
                    ll_to_data=ll_to_data, data_to_ll=data_to_ll,
                    nlines=nlines, nsamps=nsamps, fill_value=fill_value)


def determine_adjusted_data_bounds(espa_metadata, gdal_objs):
    """Determines the data boundaries specific to the NARR data resolution

    Args:
        TODO TODO TODO

    Returns:
        TODO TODO TODO
    """

    # Adjustment is specific to NARR data resolution in meters
    adjustment = 32000 * 1.0

    north_lat = float(espa_metadata.xml_object.
                      global_metadata.bounding_coordinates.north)
    south_lat = float(espa_metadata.xml_object.
                      global_metadata.bounding_coordinates.south)
    east_lon = float(espa_metadata.xml_object.
                     global_metadata.bounding_coordinates.east)
    west_lon = float(espa_metadata.xml_object.
                     global_metadata.bounding_coordinates.west)
    (ul_x, ul_y, height) = gdal_objs.ll_to_data.TransformPoint(west_lon,
                                                               north_lat)
    (ur_x, ur_y, height) = gdal_objs.ll_to_data.TransformPoint(east_lon,
                                                               north_lat)
    (lr_x, lr_y, height) = gdal_objs.ll_to_data.TransformPoint(east_lon,
                                                               south_lat)
    (ll_x, ll_y, height) = gdal_objs.ll_to_data.TransformPoint(west_lon,
                                                               south_lat)

    north_y = max(ur_y, ul_y) + adjustment
    south_y = min(ll_y, lr_y) - adjustment
    east_x = max(ur_x, lr_x) + adjustment
    west_x = min(ul_x, ll_x) - adjustment

    (ul_lon, ul_lat, height) = gdal_objs.data_to_ll.TransformPoint(west_x,
                                                                   north_y)
    (ur_lon, ur_lat, height) = gdal_objs.data_to_ll.TransformPoint(east_x,
                                                                   north_y)
    (lr_lon, lr_lat, height) = gdal_objs.data_to_ll.TransformPoint(east_x,
                                                                   south_y)
    (ll_lon, ll_lat, height) = gdal_objs.data_to_ll.TransformPoint(west_x,
                                                                   south_y)

    north_lat = max(ul_lat, ur_lat)
    south_lat = min(lr_lat, ll_lat)
    east_lon = max(ur_lon, lr_lon)
    west_lon = min(ul_lon, ll_lon)

    return DataBoundInfo(north_lat=north_lat, south_lat=south_lat,
                         east_lon=east_lon, west_lon=west_lon)


def check_within_bounds(data_bounds, longitude, latitude):
    """

    Args:
        data_bounds <DataBoundInfo>: The bounds to adhere to
        longitude <float>: The longitude to check
        latitude <float>: The latitude to check

    Returns:
        <boolean>:  True - If within the bounds
                   False - If not within the bounds
    """

    if (data_bounds.north_lat > latitude and
            data_bounds.south_lat < latitude and
            data_bounds.east_lon > longitude and
            data_bounds.west_lon < longitude):
        return True

    return False


def convert_line_to_coordinate_info(data_bounds, line):
    """Converts a line of information into a CoordinateInfo

    Args:
        data_bounds <DataBoundInfo>: The bounds to adhere to
        line <str>: The line of information read from the coordinate file

    Returns:
        <CoordinateInfo>: NamedTuple of coordinate information
    """

    (col, row, latitude, longitude) = line.split()

    ''' TODO: Should think about fixing the input file, so that this
              confusing conversion is not needed.

        When you read the file data, it is as if you are reading the values
        from the lower left to the upper right as applied to the earth.  And
        the values being read in start with an origin somewhere around the
        lower left, hence the need for the following conversion.

        NOTE: If this is changed here, then else-where in the code will break.
    '''

    # Convert to numerical types
    col = int(col)
    row = int(row)
    latitude = float(latitude)
    longitude = float(longitude)

    if longitude > 180.0:
        longitude = 360.0 - longitude
    else:
        longitude = -longitude

    if not check_within_bounds(data_bounds, longitude, latitude):
        return None
    else:
        return CoordinateInfo(col=col, row=row)


def convert_line_to_point_info(gdal_objs, min_max, line):
    """Converts a line of information into a CoordinateInfo

    Args:
        data_bounds <DataBoundInfo>: The bounds to adhere to
        line <str>: The line of information read from the coordinate file

    Returns:
        <CoordinateInfo>: NamedTuple of coordinate information
    """

    (col, row, latitude, longitude) = line.split()

    ''' TODO: Should think about fixing the input file, so that this
              confusing conversion is not needed.

        When you read the file data, it is as if you are reading the values
        from the lower left to the upper right as applied to the earth.  And
        the values being read in start with an origin somewhere around the
        lower left, hence the need for the following conversion.

        NOTE: If this is changed here, then else-where in the code will break.
    '''

    # Convert to numerical types
    col = int(col)
    row = int(row)
    latitude = float(latitude)
    longitude = float(longitude)

    if (min_max.min_row <= row and
            min_max.max_row >= row and
            min_max.min_col <= col and
            min_max.max_col >= col):

        if longitude > 180.0:
            longitude = 360.0 - longitude
        else:
            longitude = -longitude

        (map_x, map_y, height) = gdal_objs.ll_to_data.TransformPoint(longitude,
                                                                     latitude)

        return PointInfo(col=col, row=row,
                         latitude=latitude, longitude=longitude,
                         map_x=map_x, map_y=map_y)
    else:
        return None


def is_in_data(gdal_objs, point):
    """Determine if the point is located within the data
    """

    (data_x, data_y) = (util.Geo
                        .convert_mapXY_to_imageXY(point.map_x, point.map_y,
                                                  gdal_objs.data_transform))

    if data_x < 0 or data_x >= gdal_objs.nsamps:
        return False

    if data_y < 0 or data_y >= gdal_objs.nlines:
        return False

    return True


def haversine_distance(longitude_1, latitude_1, longitude_2, latitude_2):
    """Calulates the great-circle distance between two points in meters

    The points are given in decimal degrees.  Haversine formula is used.

    Args:
        longitude_1 <float>: The longitude for the first point
        latitude_1 <float>: The latitude for the first point
        longitude_2 <float>: The longitude for the second point
        latitude_2 <float>: The latitude for the second point

    Returns:
        <float>: The great-circle distance in meters

    Sources:
        https://en.wikipedia.org/wiki/Haversine_formula
    """

    # Convert to radians
    (lon_1, lat_1, lon_2, lat_2) = map(math.radians, [longitude_1, latitude_1,
                                                      longitude_2, latitude_2])

    # Figure out some sines
    sin_lon_sqrd = math.sin((lon_2 - lon_1) / 2.0)**2
    sin_lat_sqrd = math.sin((lat_2 - lat_1) / 2.0)**2

    # Compute and return the distance
    return (DIAMETER_EARTH_IN_METERS * math.asin(math.sqrt(sin_lat_sqrd +
                                                           math.cos(lat_1) *
                                                           math.cos(lat_2) *
                                                           sin_lon_sqrd)))


def grid_point_distance(gdal_objs, point, samp, line):
    """Determine the distance from the line sample to the grid point
    """

    (map_x, map_y) = (util.Geo
                      .convert_imageXY_to_mapXY(samp, line,
                                                gdal_objs.data_transform))

    (longitude, latitude, height) = gdal_objs.data_to_ll.TransformPoint(map_x,
                                                                        map_y)

    return haversine_distance(longitude, latitude,
                              point['point'].longitude,
                              point['point'].latitude)


def center_grid_point(gdal_objs, dpoints, samp, line):
    """Returns the closest point to the line sample
    """

    (map_x, map_y) = (util.Geo
                      .convert_imageXY_to_mapXY(samp, line,
                                                gdal_objs.data_transform))

    (longitude, latitude, height) = gdal_objs.data_to_ll.TransformPoint(map_x,
                                                                        map_y)

    dlist = list()
    index = 0
    for dpoint in dpoints:

        distance = haversine_distance(longitude, latitude,
                                      dpoint['point'].longitude,
                                      dpoint['point'].latitude)
        dlist.append((distance, index))
        index += 1

    # Return the closest grid point
    return sorted(dlist)[0][1]


def find_first_last_occurances(mask, lines, value):
    """Find the first and last valid data pixels for each line
    """

    for line in xrange(lines):
        l2r = mask[line]
        r2l = l2r[::-1]
        try:
            first_samp = l2r.tolist().index(value)
            last_samp = len(r2l) - 1 - r2l.tolist().index(value)
            yield (line, first_samp)
            yield (line, last_samp)
        except ValueError:
            # Some lines do not have any values and a value error is returned
            # from the index() call
            pass


def bounded_points_from_file(data_bounds, data_path):
    """Figure out the NARR points that are within the data bounds
    """

    # Generate a list of the points within the defined boundary
    with open(data_path, 'r') as coords_fd:
        points = [convert_line_to_coordinate_info(data_bounds=data_bounds,
                                                  line=line.strip())
                  for line in coords_fd]

    # Filter out the None values
    # (removes all of the way way way outside points)
    return [point for point in points if point is not None]


def gridded_points(debug, gdal_objs, lst_data_dir):
    """
    """

    # Determine full path to the file
    data_path = os.path.join(lst_data_dir, NARR_COORDINATES_FILENAME)

    # Determine buffered points
    buffered_points = bounded_points_from_file(data_bounds, data_path)

    # Determine min and max rows and cols for those within the boundary
    min_row = 9999
    max_row = -9999
    min_col = 9999
    max_col = -9999
    for point in bounded_points:
        min_row = min(min_row, point.row)
        max_row = max(max_row, point.row)
        min_col = min(min_col, point.col)
        max_col = max(max_col, point.col)

    min_max = MinMaxRowColInfo(min_row=min_row,
                               max_row=max_row,
                               min_col=min_col,
                               max_col=max_col)

    grid_rows = max_row - min_row + 1
    grid_cols = max_col - min_col + 1

    # Using the mins and maxs read the coordinates again
    # but keep a complete rectangular grid of them
    with open(data_path, 'r') as coords_fd:
        points = [convert_line_to_point_info(gdal_objs=gdal_objs,
                                             min_max=min_max,
                                             line=line.strip())
                  for line in coords_fd]

    # Filter out the None values (Only keeps the grid we want)
    points = [point for point in points if point is not None]

    if debug:
        with open('full_point_list.txt', 'w') as points_fd:
            index = 0
            for point in points:
                points_fd.write('{0},{1},{2}\n'.format(index,
                                                       point.longitude,
                                                       point.latitude))
                index += 1

    # Determine the final set of grid points
    # As well as initially mark points usable or not
    grid_points = [{'row': point.row-min_row,
                    'col': point.col-min_col,
                    'usable': is_in_data(gdal_objs=gdal_objs,
                                         point=point),
                    'point': point}
                   for point in points]

    return (grid_points, grid_rows, grid_cols)


def retrieve_narr_points(debug, gdal_objs, data_bounds, lst_data_dir):
    """Reads the NARR coordinates file and filters to a final set of points

    Args:
        data_bounds <DataBoundInfo>: The bounds to adhere to
        lst_data_dir <str>: The directory to search for the coodinate file

    Returns:
        points <list>: A list of the NARR coordinate information where each
                       element of the list is a CoordinateInfo

    Notes: The file format contains lines of the following information.
               'Grid_Column Grid_Row Grid_Latitude Grid_Longitude'
           With the following format.
               '%d %d %lf %lf'

           Each line in the file represents for the purpose of this
           application a (row, col) coodinate pair that coinsides with the
           rows and cols of the NARR data.

    """

    # Determine grid points
    (grid_points, grid_rows, grid_cols) = gridded_points(debug,
                                                         gdal_objs,
                                                         lst_data_dir)

    raster_data = (gdal_objs.data_ds.GetRasterBand(1)
                   .ReadAsArray(0, 0, gdal_objs.nsamps, gdal_objs.nlines))

    mask = np.zeros(shape=raster_data.shape, dtype=np.bool)

    # Set all the valid data points to True
    mask[raster_data != gdal_objs.fill_value] = True
    del raster_data

    # Process through the mask and generate pairs for the east/west edges
    ew_edges = sorted([pair
                       for pair
                       in find_first_last_occurances(mask, gdal_objs.nlines,
                                                     True)])

    # Add all the data points between the first and last line edges
    first_line = ew_edges[:2]
    last_line = ew_edges[-2:]
    min_line = first_line

    ew_edges.extend([(first_line[0][0], samp)
                     for samp in xrange(first_line[0][1] + 1,
                                        first_line[1][1])])

    ew_edges.extend([(last_line[0][0], samp)
                     for samp in xrange(last_line[0][1] + 1,
                                        last_line[1][1])])

    # Remove duplicates and re-sort
    ew_edges = sorted(list(set(ew_edges)))

    for pair in ew_edges:
        cc_point = center_grid_point(gdal_objs, grid_points, pair[1], pair[0])

        '''
            UL UC UR
            CL CC CR
            LL LC LR
        '''

        # Figure out all the remaining possible grid points
        ul_point = cc_point + grid_cols - 1
        uc_point = ul_point + 1
        ur_point = uc_point + 1
        cl_point = cc_point - 1
        cr_point = cc_point + 1
        ll_point = cc_point - grid_cols - 1
        lc_point = ll_point + 1
        lr_point = lc_point + 1

        # Get all of the distance to the outer grid points
        ul_dist = grid_point_distance(gdal_objs, grid_points[ul_point],
                                      pair[1], pair[0])
        uc_dist = grid_point_distance(gdal_objs, grid_points[uc_point],
                                      pair[1], pair[0])
        ur_dist = grid_point_distance(gdal_objs, grid_points[ur_point],
                                      pair[1], pair[0])

        cl_dist = grid_point_distance(gdal_objs, grid_points[cl_point],
                                      pair[1], pair[0])
        cr_dist = grid_point_distance(gdal_objs, grid_points[cr_point],
                                      pair[1], pair[0])

        ll_dist = grid_point_distance(gdal_objs, grid_points[ll_point],
                                      pair[1], pair[0])
        lc_dist = grid_point_distance(gdal_objs, grid_points[lc_point],
                                      pair[1], pair[0])
        lr_dist = grid_point_distance(gdal_objs, grid_points[lr_point],
                                      pair[1], pair[0])

        # Determine quadrant by using the average quadrant distances

        avg_dist_ul = (cl_dist + ul_dist + uc_dist) / 3.0
        avg_dist_ur = (uc_dist + ur_dist + cr_dist) / 3.0
        avg_dist_lr = (cr_dist + lr_dist + lc_dist) / 3.0
        avg_dist_ll = (lc_dist + ll_dist + cl_dist) / 3.0

        min_dist = min(avg_dist_ul, avg_dist_ur, avg_dist_lr, avg_dist_ll)

        if min_dist == avg_dist_ul:
            grid_points[cl_point]['usable'] = True
            grid_points[ul_point]['usable'] = True
            grid_points[uc_point]['usable'] = True
        elif min_dist == avg_dist_ur:
            grid_points[uc_point]['usable'] = True
            grid_points[ur_point]['usable'] = True
            grid_points[cr_point]['usable'] = True
        elif min_dist == avg_dist_lr:
            grid_points[cr_point]['usable'] = True
            grid_points[lr_point]['usable'] = True
            grid_points[lc_point]['usable'] = True
        else:
            grid_points[lc_point]['usable'] = True
            grid_points[ll_point]['usable'] = True
            grid_points[cl_point]['usable'] = True

    return grid_points


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

    # Figure out the gdal objects
    gdal_objs = initialize_gdal_objects(espa_metadata=espa_metadata)

    # Determin NARR adjusted data bounds
    data_bounds = determine_adjusted_data_bounds(espa_metadata=espa_metadata,
                                                 gdal_objs=gdal_objs)
    logger.debug(str(data_bounds))

    # Coordinate Information
    points = retrieve_narr_points(debug=args.debug,
                                  gdal_objs=gdal_objs,
                                  data_bounds=data_bounds,
                                  lst_data_dir=args.lst_data_dir)

    if args.debug:
        with open('usable_point_list.txt', 'w') as usable_fd:
            with open('unusable_point_list.txt', 'w') as unusable_fd:
                index = 0
                for point in points:
                    line = ('{0},{1},{2},{3},{4}\n'
                            .format(index, point['row'], point['col'],
                                    point['point'].longitude,
                                    point['point'].latitude))

                    if point['usable']:
                        usable_fd.write(line)
                    if not point['usable']:
                        unusable_fd.write(line)

                index += 1

    del points


if __name__ == '__main__':
    main()
