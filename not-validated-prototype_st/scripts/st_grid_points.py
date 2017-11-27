
'''
    File: st_build_points.py

    Purpose: Builds a directory structure of points and required information
             to be used for generation and interpolation of atmospheric
             information by follow on applications.

    Project: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    License: NASA Open Source Agreement 1.3
'''

import logging
import array
import struct
from collections import namedtuple


# Grid Point Information
PointInfo = namedtuple('PointInfo',
                       ('col', 'row', 'lat', 'lon', 'map_y', 'map_x'))

GridPointInfo = namedtuple('GridPointInfo',
                           ('index',
                            'run_modtran',
                            'col', 'row',
                            'narr_col', 'narr_row',
                            'lat', 'lon',
                            'map_y', 'map_x'))

GRID_POINT_FMT = 'HBBBHHffff'
GRID_POINT_HEADER_NAME = 'grid_points.hdr'
GRID_POINT_BINARY_NAME = 'grid_points.bin'


def write_grid_points(grid_points, grid_rows, grid_cols):
    """Writes grid points to a binary file with a header

    Binary is used so tha precision in the floats is not lost.

    Args:
        grid_points [<PointInfo>]: The list of grid points
    """

    with open(GRID_POINT_HEADER_NAME, 'w') as ascii_fd:
        ascii_fd.write('{}\n'.format(len(grid_points)))
        ascii_fd.write('{}\n'.format(grid_rows))
        ascii_fd.write('{}\n'.format(grid_cols))

    with open(GRID_POINT_BINARY_NAME, 'wb') as binary_fd:
        point_struct = struct.Struct(GRID_POINT_FMT)

        # Allocate enough buffer space for all the points we will write
        buffer_size = '\0' * point_struct.size * len(grid_points)
        point_buffer = array.array('c', buffer_size)

        index = 0
        for point in grid_points:
            point_struct.pack_into(point_buffer, index * point_struct.size,
                                   point['index'],
                                   int(point['run_modtran']),
                                   point['row'],
                                   point['col'],
                                   point['narr_row'],
                                   point['narr_col'],
                                   point['point'].lon,
                                   point['point'].lat,
                                   point['point'].map_x,
                                   point['point'].map_y)
            index += 1

        binary_fd.write(point_buffer)


def read_grid_points():
    """Read grid points from the binary file back into a structure

    Returns:
        grid_points [<GridPointInfo>]: The list of grid points
    """

    logger = logging.getLogger(__name__)

    count = 0
    grid_rows = 0
    grid_cols = 0
    with open(GRID_POINT_HEADER_NAME, 'r') as ascii_fd:
        count = int(ascii_fd.readline().strip())
        grid_rows = int(ascii_fd.readline().strip())
        grid_cols = int(ascii_fd.readline().strip())

    logger.info('Reading [{}] points from grid file'.format(count))

    grid_points = list()
    with open(GRID_POINT_BINARY_NAME, 'rb') as binary_fd:
        point_struct = struct.Struct(GRID_POINT_FMT)

        for position in xrange(count):
            logger.debug('Reading point {}'.format(position))
            point_buffer = bytearray(binary_fd.read(point_struct.size))
            (index, run_modtran, row, col, narr_row, narr_col,
             lon, lat, map_x, map_y) = point_struct.unpack(point_buffer)

            grid_points.append(GridPointInfo(index=index,
                                             run_modtran=run_modtran,
                                             row=row, col=col,
                                             narr_row=narr_row,
                                             narr_col=narr_col,
                                             lon=lon, lat=lat,
                                             map_x=map_x, map_y=map_y))

    return (grid_points, grid_rows, grid_cols)
