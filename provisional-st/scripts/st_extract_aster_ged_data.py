#! /usr/bin/env python3

'''
    FILE: st_extract_aster_ged_data.py

    PURPOSE: Extract the data used by Surface Temperature from ASTER GED 
             tiles.  The results are meant to be used for a trimmed local 
             copy of the ASTER GED product.

             This is an initial draft of the tool used to assess different
             options.  It should ultimately be moved to the auxiliary tools.

    PROJECT: Land Satellites Data Systems (LSDS) Science Research and
             Development (LSRD) at the USGS EROS

    LICENSE: NASA Open Source Agreement 1.3
'''

import os
import sys
import logging


import numpy as np
from osgeo import gdal, osr


from st_exceptions import InaccessibleTileError


# Import local modules
import st_utilities as util
import emissivity_utilities as emis_util


def extract_aster_data(url, filename, intermediate):
    """Extracts the internal band(s) data for later processing

    Args:
        url <str>: URL to retrieve the file from
        filename <str>: Base HDF filename to extract from
        intermediate <bool>: Keep any intermediate products generated

    Returns:
        <numpy.2darray>: SDev Band 13 data
        <numpy.2darray>: SDev Band 14 data
        <int>: Samples in the data
        <int>: Lines in the data
        <2x3:float>: GDAL Affine transformation matrix
                     [0] - Map X of upper left corner
                     [1] - Pixel size in X direction
                     [2] - Y rotation
                     [3] - Map Y of upper left corner
                     [4] - X rotation
                     [5] - Pixel size in Y direction
        <bool>: True if the ASTER tile is available, False otherwise
    """

    logger = logging.getLogger(__name__)

    # Build the HDF5 filename for the tile
    h5_file_path = ''.join([filename, '.h5'])

    if not os.path.exists(h5_file_path):
        emis_util.download_aster_ged_tile(url=url, h5_file_path=h5_file_path)

    # There are cases where the emissivity data will not be available
    # (for example, in water regions).
    aster_b13_mean_data = []
    aster_b14_mean_data = []
    aster_b13_sdev_data = []
    aster_b14_sdev_data = []
    aster_ndvi_mean_data = []
    samps = 0
    lines = 0
    geo_transform = []
    if not os.path.exists(h5_file_path):
        # The ASTER tile is not available, so don't try to process it
        return (aster_b13_mean_data, aster_b14_mean_data, aster_b13_sdev_data, 
                aster_b14_sdev_data, aster_ndvi_mean_data, aster_lat_data, 
                aster_lon_data, samps, lines, geo_transform, False)

    # Define the sub-dataset names
    emis_mean_ds_name = ''.join(['HDF5:"', h5_file_path, '"://Emissivity/Mean'])
    emis_sdev_ds_name = ''.join(['HDF5:"', h5_file_path, '"://Emissivity/SDev'])
    ndvi_mean_ds_name = ''.join(['HDF5:"', h5_file_path, '"://NDVI/Mean'])
    lat_ds_name = ''.join(['HDF5:"', h5_file_path, '"://Geolocation/Latitude'])
    lon_ds_name = ''.join(['HDF5:"', h5_file_path, '"://Geolocation/Longitude'])

    logger.debug(lat_ds_name)
    logger.debug(lon_ds_name)

    aster_b13_mean_data = util.Dataset.extract_raster_data(emis_mean_ds_name, 4)
    aster_b14_mean_data = util.Dataset.extract_raster_data(emis_mean_ds_name, 5)
    aster_b13_sdev_data = util.Dataset.extract_raster_data(emis_sdev_ds_name, 4)
    aster_b14_sdev_data = util.Dataset.extract_raster_data(emis_sdev_ds_name, 5)
    aster_ndvi_mean_data = util.Dataset.extract_raster_data(ndvi_mean_ds_name,1)
    aster_lat_data = util.Dataset.extract_raster_data(lat_ds_name, 1)
    aster_lon_data = util.Dataset.extract_raster_data(lon_ds_name, 1)

    # Determine the minimum and maximum latitude and longitude
    x_min = aster_lon_data.min()
    x_max = aster_lon_data.max()
    y_min = aster_lat_data.min()
    y_max = aster_lat_data.max()

    # Determine the resolution and dimensions of the ASTER data
    (x_res, y_res, samps, lines) = (
        emis_util.data_resolution_and_size(lat_ds_name,
                                           x_min, x_max, y_min, y_max))

    # Remove the HDF5 tile since we no longer need it
    if not intermediate:
        if os.path.exists(h5_file_path):
            os.unlink(h5_file_path)

    # Build the geo transform
    geo_transform = [x_min, x_res, 0, y_max, 0, -y_res]

    return (aster_b13_mean_data, aster_b14_mean_data, aster_b13_sdev_data, 
            aster_b14_sdev_data, aster_ndvi_mean_data, aster_lat_data, 
            aster_lon_data, samps, lines, geo_transform, True)


def generate_emis_tiles(tile_name_b13, tile_name_b14, aster_b13_data, 
                        aster_b14_data, samps, lines, transform, wkt, 
                        no_data_value):
    """Generate b13 and b14 emissivity tiles

    Args:
        tile_name_b13 <str>: Filename to create for the b13 tile
        tile_name_b14 <str>: Filename to create for the b14 tile
        aster_b13_data <numpy.2darray>: Unscaled band 13 ASTER tile data
        aster_b14_data <numpy.2darray>: Unscaled band 14 ASTER tile data
        samps <int>: Samples in the data
        lines <int>: Lines in the data
        transform <2x3:float>: GDAL Affine transformation matrix
                               [0] - Map X of upper left corner
                               [1] - Pixel size in X direction
                               [2] - Y rotation
                               [3] - Map Y of upper left corner
                               [4] - X rotation
                               [5] - Pixel size in Y direction
        wkt <str>: Well-Known-Text describing the projection
        no_data_value <float>: Value to use for fill
    """

    logger = logging.getLogger(__name__)

    # Create the estimated Landsat EMIS b13 raster output tile
    logger.info('Creating an estimated Landsat b13 EMIS '
                'tile {}'.format(tile_name_b13))
    util.Geo.generate_raster_file(gdal.GetDriverByName('GTiff'),
                                  tile_name_b13,
                                  aster_b13_data,
                                  samps, lines,
                                  transform,
                                  wkt,
                                  no_data_value,
                                  gdal.GDT_Int16)

    # Create the estimated Landsat EMIS b14 raster output tile
    logger.info('Creating an estimated Landsat b14 EMIS '
                'tile {}'.format(tile_name_b14))
    util.Geo.generate_raster_file(gdal.GetDriverByName('GTiff'),
                                  tile_name_b14,
                                  aster_b14_data,
                                  samps, lines,
                                  transform,
                                  wkt,
                                  no_data_value,
                                  gdal.GDT_Int16)


def generate_aster_ndvi_tile(tile_name, ndvi_data, samps, lines, transform,
                             wkt, no_data_value):
    """Generate ASTER NDVI tile

    Args:
        tile_name <str>: Filename to create for the tile
        ndvi_data <numpy.2darray>: Unscaled NDVI Band ASTER tile data
        samps <int>: Samples in the data
        lines <int>: Lines in the data
        transform <2x3:float>: GDAL Affine transformation matrix
                               [0] - Map X of upper left corner
                               [1] - Pixel size in X direction
                               [2] - Y rotation
                               [3] - Map Y of upper left corner
                               [4] - X rotation
                               [5] - Pixel size in Y direction
        wkt <str>: Well-Known-Text describing the projection
        no_data_value <float>: Value to use for fill
    """

    logger = logging.getLogger(__name__)

    # Create the ASTER NDVI raster output tile
    logger.info('Creating an ASTER NDVI tile {}'.format(tile_name))
    util.Geo.generate_raster_file(gdal.GetDriverByName('GTiff'),
                                  tile_name,
                                  ndvi_data,
                                  samps, lines,
                                  transform,
                                  wkt,
                                  no_data_value,
                                  gdal.GDT_Int16)


def generate_emis_stdev_tiles(tile_name_b13, tile_name_b14, 
                             aster_b13_stdev_data, aster_b14_stdev_data, 
                             samps, lines, transform, wkt, no_data_value):
    """Generate ASTER emissivity standard deviation tiles

    Args:
        tile_name_b13 <str>: Filename to create for the b13 tile
        tile_name_b14 <str>: Filename to create for the b14 tile
        aster_b13_data <numpy.2darray>: Unscaled band 13 ASTER stdev tile data
        aster_b14_data <numpy.2darray>: Unscaled band 14 ASTER stdev tile data
        samps <int>: Samples in the data
        lines <int>: Lines in the data
        transform <2x3:float>: GDAL Affine transformation matrix
                               [0] - Map X of upper left corner
                               [1] - Pixel size in X direction
                               [2] - Y rotation
                               [3] - Map Y of upper left corner
                               [4] - X rotation
                               [5] - Pixel size in Y direction
        wkt <str>: Well-Known-Text describing the projection
        no_data_value <float>: Value to use for fill
    """

    logger = logging.getLogger(__name__)

    # Create the ASTER EMIS b13 standard deviation raster output tile
    logger.info('Creating ASTER EMIS b13 STDEV tile {}'.format(tile_name_b13))
    util.Geo.generate_raster_file(gdal.GetDriverByName('GTiff'),
                                  tile_name_b13,
                                  aster_b13_stdev_data,
                                  samps, lines,
                                  transform,
                                  wkt,
                                  no_data_value,
                                  gdal.GDT_Int16)

    # Create the ASTER EMIS b14 standard deviation raster output tile
    logger.info('Creating ASTER EMIS b14 STDEV tile {}'.format(tile_name_b14))
    util.Geo.generate_raster_file(gdal.GetDriverByName('GTiff'),
                                  tile_name_b14,
                                  aster_b14_stdev_data,
                                  samps, lines,
                                  transform,
                                  wkt,
                                  no_data_value,
                                  gdal.GDT_Int16)


def generate_tiles(st_data_dir, url, wkt, no_data_value, intermediate):
    """Generate tiles for emissivity standard deviation from ASTER data

    Args:
        st_data_dir <str>: Location of the ST data files
        url <str>: URL to retrieve the file from
        wkt <str>: Well-Known-Text describing the projection
        no_data_value <float>: Value to use for fill
        intermediate <bool>: Keep any intermediate products generated

    Returns:
        list(<str>): Standard deviation emissivity tile names
    """

    logger = logging.getLogger(__name__)

    # Read the ASTER GED tile list
    ged_tile_file = 'aster_ged_tile_list.txt'
    with open(os.path.join(st_data_dir, ged_tile_file)) as ged_file:
        tiles = [os.path.splitext(line.rstrip('\n'))[0] for line in ged_file]

    for filename in tiles:

        # Build the output tile names
        ls_emis_mean_tile_name_b13 = ''.join([filename, '_emis_b13.tif'])
        ls_emis_mean_tile_name_b14 = ''.join([filename, '_emis_b14.tif'])
        ls_emis_stdev_tile_name_b13 = ''.join([filename, '_emis_stdev_b13.tif'])
        ls_emis_stdev_tile_name_b14 = ''.join([filename, '_emis_stdev_b14.tif'])
        aster_ndvi_mean_tile_name = ''.join([filename, '_ndvi.tif'])

        # Read the ASTER data
        (aster_b13_mean_data, aster_b14_mean_data, aster_b13_stdev_data, 
         aster_b14_stdev_data, aster_ndvi_mean_data, aster_lat_data,
         aster_lon_data, samps, lines, transform, aster_data_available) = (
             extract_aster_data(url=url,
                                filename=filename,
                                intermediate=intermediate))

        # Fail if a tile can't be read, but it is in the ASTER GED
        if not aster_data_available:
            raise InaccessibleTileError(
                'Cannot reach tile {} in ASTER GED'.format(filename))

        # Build the emissivity mean tiles
        generate_emis_tiles(tile_name_b13=ls_emis_mean_tile_name_b13,
                            tile_name_b14=ls_emis_mean_tile_name_b14,
                            aster_b13_data=aster_b13_mean_data,
                            aster_b14_data=aster_b14_mean_data,
                            samps=samps,
                            lines=lines,
                            transform=transform,
                            wkt=wkt,
                            no_data_value=no_data_value)

        del aster_b13_mean_data
        del aster_b14_mean_data

        # Build the NDVI tile
        generate_aster_ndvi_tile(tile_name=aster_ndvi_mean_tile_name,
                                 ndvi_data=aster_ndvi_mean_data,
                                 samps=samps,
                                 lines=lines,
                                 transform=transform,
                                 wkt=wkt,
                                 no_data_value=no_data_value)

        del aster_ndvi_mean_data

        # Build the emissivity standard deviation tiles
        generate_emis_stdev_tiles(tile_name_b13=ls_emis_stdev_tile_name_b13,
                                  tile_name_b14=ls_emis_stdev_tile_name_b14,
                                  aster_b13_stdev_data=aster_b13_stdev_data,
                                  aster_b14_stdev_data=aster_b14_stdev_data,
                                  samps=samps,
                                  lines=lines,
                                  transform=transform,
                                  wkt=wkt,
                                  no_data_value=no_data_value)

        del aster_b13_stdev_data
        del aster_b14_stdev_data


def build_ls_emis_data(server_name, server_path, st_data_dir, 
                       no_data_value, intermediate):
    """Build estimated Landsat Emissivity Data

    Args:
        server_name <str>: Name of the ASTER GED server
        server_path <str>: Path on the ASTER GED server
        st_data_dir <str>: Location of the ST data files
        no_data_value <int>: No data (fill) value to use
        intermediate <bool>: Keep any intermediate products generated
    """

    logger = logging.getLogger(__name__)

    # Specify the base URL to use for retrieving the ASTER GED data
    url = ''.join(['http://', server_name, server_path])

    # The ASTER data is in geographic projection so specify that here
    ds_srs = osr.SpatialReference()
    ds_srs.ImportFromEPSG(4326)
    geographic_wkt = ds_srs.ExportToWkt()

    generate_tiles(st_data_dir=st_data_dir,
                   url=url,
                   wkt=geographic_wkt,
                   no_data_value=no_data_value,
                   intermediate=intermediate)



def main():
    """Extract ASTER GED tiles, pull out the bands that ST needs, and convert 
       them to GeoTIFF.
    """

    args = emis_util.retrieve_command_line_arguments()

    # Check logging level
    debug_level = logging.INFO

    if args.debug:
        debug_level = logging.DEBUG

    # Configure logging
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:%(funcName)s'
                                ' -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=debug_level)

    logger = logging.getLogger(__name__)

    logger.info('*** Extracting ASTER GED Tiles ***')

    try:
        # Register all the gdal drivers
        gdal.AllRegister()

        # Get the data directory from the environment
        st_data_dir = emis_util.get_env_var('ST_DATA_DIR', None)

        # Call the main processing routine
        build_ls_emis_data(server_name=args.aster_ged_server_name,
                           server_path=args.aster_ged_server_path,
                           st_data_dir=st_data_dir,
                           no_data_value=util.INTERMEDIATE_NO_DATA_VALUE,
                           intermediate=args.intermediate)
    except Exception:
        logger.exception('Processing failed')
        sys.exit(1)  # EXIT FAILURE

    logger.info('*** Extract ASTER GED Tiles - Complete ***')


if __name__ == '__main__':
    main()
