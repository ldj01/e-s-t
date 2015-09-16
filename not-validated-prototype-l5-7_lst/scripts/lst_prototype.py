#! /usr/bin/env python


import os
import sys
import logging
import copy
import math
from argparse import ArgumentParser


import metadata_api


# ============================================================================
class NARR_Point(object):
    '''
    Description:
        An object structure to hold associated values.
    '''

    row = None        # Ones based indexing from the source file
    col = None        # Ones based indexing from the source file
    latitude = None
    longitude = None
    northing = None
    easting = None


# ============================================================================
class Points_Info(object):
    '''
    Description:
        An object structure to hold associated values.
    '''

    min_row = None
    max_row = None
    min_col = None
    max_col = None
    num_rows = None
    num_cols = None
    num_points = None
    grid_points = None


# ============================================================================
def read_narr_coordinates(base_data_dir):
    '''
    Description:
        Read the NARR coordinates into memory.

    Returns:
        grid_data - The list of NARR points with longitude fixed.
    '''

    narr_path = os.path.join(base_data_dir, 'narr_coordinates.txt')

    grid_data = list()
    with open(narr_path, 'r') as narr_fd:
        for line in narr_fd:
            line = line.strip()
            (grid_col, grid_row, grid_lat, grid_lon) = line.split()

            grid_lon = float(grid_lon)

            point = NARR_Point()
            point.row = int(grid_row)
            point.col = int(grid_col)
            point.latitude = float(grid_lat)

            #  When you read the file data, it is as if you are reading the
            #  values from the lower left to the upper right as applied to
            #  the earth.  And the values being read in start with an origin
            #  somewhere around the lower left, hence the need for the
            #  following conversion.
            #
            #  NOTE - If this is changed here, then else-where in the code
            #         will break.
            if grid_lon > 180.0:
                point.longitude = 360.0 - grid_lon
            else:
                point.longitude = -grid_lon

            grid_data.append(point)

    return grid_data


# ============================================================================
def convert_latitude_longitude_to_utm_northing_easting(gm, grid_points):
    '''
    Description:
        TODO TODO TODO
    '''

"""
void convert_ll_to_utm
(
    Input_t *input,
    REANALYSIS_POINTS *points /* I/O: The coordinate points to be used */
)
{
    int point;
    double a = UTM_EQUATORIAL_RADIUS; /* equatorial radius */
    double b = UTM_POLAR_RADIUS;      /* polar radius */
    double k0 = UTM_SCALE_FACTOR;     /* scale factor */
    double ecc;                       /* eccentricity */
    double ecc_prime_sqrd;            /* prime of eccentricity squared */
    double n;
    double nu;
    double a0;
    double b0;
    double c0;
    double d0;
    double e0;
    double ki;
    double kii;
    double kiii;
    double kiv;
    double kv;
    double zone_cm;
    double p;
    double lat_rad;
    double s;

    /* Calculate zone central meridian in degrees and radians */
    zone_cm = (6.0 * input->meta.zone) - 183.0;

    ecc = sqrt (1.0 - ((b / a) * (b / a)));
    ecc_prime_sqrd = (ecc * ecc) / (1.0 - ecc * ecc);
    n = (a - b) / (a + b);

    /* Calculate meridional arc length */
    a0 = a * (1.0 - n
              + (5.0 * n * n / 4.0) * (1.0 - n)
              + (81.0 * pow (n, 4) / 64.0) * (1.0 - n));

    b0 = (3.0 * a * n / 2.0) * (1.0 - n
                                - (7.0 * n * n / 8.0) * (1.0 - n)
                                + 55.0 * pow (n, 4) / 64.0);

    c0 = (15.0 * a * n * n / 16.0) * (1.0 - n
                                      + (3.0 * n * n / 4.0) * (1.0 - n));

    d0 = (35 * a * pow (n, 3) / 48.0) * (1.0 - n
                                         + 11.0 * n * n / 16.0);

    e0 = (315.0 * a * pow (n, 4) / 51.0) * (1.0 - n);

    for (point = 0; point < points->num_points; point++)
    {
        /* Distance from the point to the zone central meridian in radians */
        p = (points->lon[point] - zone_cm) * RADIANS_PER_DEGREE;

        /* Convert latitude for the point to radians */
        lat_rad = points->lat[point] * RADIANS_PER_DEGREE;

        nu = a / sqrt (1.0 - (ecc * sin (lat_rad))
                             * (ecc * sin (lat_rad)));

        s = a0 * lat_rad
            - b0 * sin (2 * lat_rad)
            + c0 * sin (4 * lat_rad)
            - d0 * sin (6 * lat_rad)
            + e0 * sin (8 * lat_rad);

        /* Coefficients for UTM coordinates */
        ki = s * k0;

        kii = nu * sin (lat_rad) * cos (lat_rad) * k0 / 2.0;

        kiii = (nu * sin (lat_rad) * pow (cos (lat_rad), 3) / 24.0)
               * (5.0 - tan (lat_rad) * tan (lat_rad)
                  + 9.0 * ecc_prime_sqrd * pow (cos (lat_rad), 2)
                  + 4.0 * ecc_prime_sqrd * ecc_prime_sqrd
                    * pow (cos (lat_rad), 4))
               * k0;

        kiv = nu * cos (lat_rad) * k0;

        kv = pow (cos (lat_rad), 3)
             * (nu / 6.0)
             * (1.0 - tan (lat_rad) * tan (lat_rad)
                + ecc_prime_sqrd * cos (lat_rad) * cos (lat_rad))
             * k0;

        /* Calculate UTM coordinates */
        points->utm_easting[point] = UTM_FALSE_EASTING
                                     + (kiv * p + kv * pow (p, 3));
        points->utm_northing[point] = (ki + kii * p * p + kiii * pow (p, 4));
    }
}
"""


# ============================================================================
def build_points(gm, base_data_dir):
    '''
    Description:
        TODO TODO TODO
    '''

    logger = logging.getLogger(__name__)

    points = Points_Info()

    '''
    expand range to include NARR points outside image for edge pixels */
    TODO - 0.2deg at higher latitudes will not be sufficient for the
           longitudinal test, since at lat(72deg) 1deg lon = 34504.22meters
           and the NARR data is 34k between points.

           This is probably only a CONUS quick and dirty solution.

    NOTE - MERRA is even farther apart so this will not work for that. */
    '''
    buffered_north_lat = float(gm.bounding_coordinates.north) + 0.2
    buffered_south_lat = float(gm.bounding_coordinates.south) - 0.2
    buffered_east_lat = float(gm.bounding_coordinates.east) + 0.2
    buffered_west_lat = float(gm.bounding_coordinates.west) - 0.2
    logger.info('Buffered Latitude and Longitude for point determination')
    logger.info('  buffered north latitude = {0}'.format(buffered_north_lat))
    logger.info('  buffered south latitude = {0}'.format(buffered_south_lat))
    logger.info('  buffered east longitude = {0}'.format(buffered_east_lat))
    logger.info('  buffered west longitude = {0}'.format(buffered_west_lat))

    try:
        grid_data = read_narr_coordinates(base_data_dir)
    except Exception:
        logger.error('Failed reading auxillary NARR dataset')
        raise

    '''
    Determine the points in the NARR dataset that fall within our buffered
    Landsat area using logical operators lessThanLat and greaterThanLat are
    values where the NARR values are less than or greater than the edges of
    the Landsat corners values respectively pixels that are true in both fall
    within the Landsat scene the same thing is done with longitude values
    '''
    min_row = 1000
    max_row = -1000
    min_col = 1000
    max_col = -1000
    for point in grid_data:
        if (buffered_north_lat > point.latitude and
                buffered_south_lat < point.latitude and
                buffered_west_lat < point.longitude and
                buffered_east_lat > point.longitude):
            min_row = min(min_row, point.row)
            max_row = max(max_row, point.row)
            min_col = min(min_col, point.col)
            max_col = max(max_col, point.col)

    # The row/col values from the file are ones based, and we need zero based,
    # so fix them
    points.min_row = min_row - 1
    points.max_row = max_row - 1
    points.min_col = min_col - 1
    points.max_col = max_col - 1
    points.num_rows = max_row - min_row + 1
    points.num_cols = max_col - min_col + 1
    logger.info('NARR Points (Min/Max) Range')
    logger.info('  min_row = {0}'.format(points.min_row))
    logger.info('  max_row = {0}'.format(points.max_row))
    logger.info('  min_col = {0}'.format(points.min_col))
    logger.info('  max_col = {0}'.format(points.max_col))
    logger.info('  num_rows = {0}'.format(points.num_rows))
    logger.info('  num_cols = {0}'.format(points.num_cols))

    points.grid_points = list()
    # Retain only the points within the rectangle
    for row in xrange(points.min_row, points.max_row + 1):
        for col in xrange(points.min_col, points.max_col + 1):
            index = (row - min_row) * points.num_cols + (col - min_col)

            point = copy.deepcopy(grid_data[index])
            points.grid_points.append(point)
    del (grid_data)

    # TODO TODO TODO
    # /* Convert lat/lon to UTM northing/easting */
    # convert_latitude_longitude_to_utm_northing_easting(gm,
    # TODO TODO TODO

    return points


# ============================================================================
def process(xml_filename, base_data_dir):
    '''
    Description:
        TODO TODO TODO
    '''

    logger = logging.getLogger(__name__)

    # Read the XML metadata
    espa_xml = metadata_api.parse(xml_filename, silence=True)
    # Grab the global metadata object
    gm = espa_xml.get_global_metadata()

    # Build the points to use from the auxillary NARR dataset
    points = build_points(gm, base_data_dir)

    # TODO TODO TODO
    # TODO TODO TODO
    # TODO TODO TODO
    # TODO TODO TODO
    # TODO TODO TODO
    # TODO TODO TODO

    del (points)


# ============================================================================
if __name__ == '__main__':
    '''
    Description:
        Read the NARR coordinates into memory.
    '''

    # Build the command line argument parser
    description = ('Retrieve ASTER data application')
    parser = ArgumentParser(description=description)

    # ---- Add parameters ----
    # Required parameters
    parser.add_argument('--xml',
                        action='store', dest='xml_filename', required=True,
                        help='The XML metadata file to use')

    # Parse the command line arguments
    args = parser.parse_args()

    # Configure logging
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:%(funcName)s'
                                ' -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)

    logger = logging.getLogger(__name__)

    base_data_dir = os.environ.get('LST_DATA_DIR')
    if base_data_dir is None:
        logger.info('Missing environment variable LST_DATA_DIR')
        sys.exit(1)  # EXIT FAILURE

    if not os.path.isdir(base_data_dir):
        logger.info('LST_DATA_DIR directory does not exist')
        sys.exit(1)  # EXIT FAILURE

    try:
        # Call the main processing routine
        process(args.xml_filename, base_data_dir)

    except Exception:
        logger.exception('Processing failed')
        sys.exit(1)  # EXIT FAILURE

    sys.exit(0)  # EXIT SUCCESS
