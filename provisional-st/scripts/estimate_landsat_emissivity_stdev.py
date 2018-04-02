#! /usr/bin/env python

'''
    FILE: estimate_landsat_emissivity_stdev.py

    PURPOSE: Estimates a Landsat Emissivity standard deviation product from
             ASTER Emissivity. The results are meant to be used for generation
             of a Surface Temperature product.

    PROJECT: Land Satellites Data Systems (LSDS) Science Research and
             Development (LSRD) at the USGS EROS

    LICENSE: NASA Open Source Agreement 1.3

    Algorithm Authors:

        Glynn C. Hulley
        Research Scientist
        NASA Jet Propulsion Laboratory
        email: glynn.hulley@nasa.gov

        Simon J. Hook
        Research Scientist
        NASA Jet Propulsion Laboratory
        email: Simon.J.Hook@jpl.nasa.gov
'''

import os
import sys
import logging


import numpy as np
from osgeo import gdal, osr


from espa import Metadata
from st_exceptions import NoTilesError, InaccessibleTileError


# Import local modules
import st_utilities as util
import emissivity_utilities as emis_util


ASTER_GED_N_FORMAT = 'AG100.v003.{0:02}.{1:04}.0001'
ASTER_GED_P_FORMAT = 'AG100.v003.{0:02}.{1:03}.0001'

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
    aster_b13_sdev_data = []
    aster_b14_sdev_data = []
    samps = 0
    lines = 0
    geo_transform = []
    if not os.path.exists(h5_file_path):
        # The ASTER tile is not available, so don't try to process it
        return (aster_b13_sdev_data, aster_b14_sdev_data, samps, lines,
                geo_transform, False)

    # Define the sub-dataset names
    emis_sdev_ds_name = ''.join(['HDF5:"', h5_file_path,
                                 '"://Emissivity/SDev'])
    lat_ds_name = ''.join(['HDF5:"', h5_file_path,
                           '"://Geolocation/Latitude'])
    lon_ds_name = ''.join(['HDF5:"', h5_file_path,
                           '"://Geolocation/Longitude'])

    logger.debug(lat_ds_name)
    logger.debug(lon_ds_name)

    aster_b13_sdev_data = emis_util.extract_raster_data(emis_sdev_ds_name, 4)
    aster_b14_sdev_data = emis_util.extract_raster_data(emis_sdev_ds_name, 5)
    aster_lat_data = emis_util.extract_raster_data(lat_ds_name, 1)
    aster_lon_data = emis_util.extract_raster_data(lon_ds_name, 1)

    # Determine the minimum and maximum latitude and longitude
    x_min = aster_lon_data.min()
    x_max = aster_lon_data.max()
    y_min = aster_lat_data.min()
    y_max = aster_lat_data.max()

    del aster_lon_data
    del aster_lat_data

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

    return (aster_b13_sdev_data, aster_b14_sdev_data, samps, lines,
            geo_transform, True)


def generate_emis_stdev_tile(tile_name, aster_b13_stdev_data,
                             aster_b14_stdev_data, samps, lines, transform,
                             wkt, no_data_value):
    """Generate ASTER emissivity standard deviation values for the tile

    Args:
        tile_name <str>: Filename to create for the tile
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

    # Save the no data locations.
    aster_b13_stdev_no_data_locations = np.where(aster_b13_stdev_data ==
                                                 no_data_value)

    # Scale the data
    aster_b13_stdev_data = aster_b13_stdev_data * 0.0001

    # Save the no data locations.
    aster_b14_stdev_no_data_locations = np.where(aster_b14_stdev_data ==
                                                 no_data_value)

    # Scale the data
    aster_b14_stdev_data = aster_b14_stdev_data * 0.0001

    # Create the estimated Landsat EMIS stdev data.
    emis_stdev_data = np.sqrt((aster_b13_stdev_data**2
                               + aster_b14_stdev_data**2)/2)

    # Re-apply the no data locations.
    emis_stdev_data[aster_b13_stdev_no_data_locations] = no_data_value
    emis_stdev_data[aster_b14_stdev_no_data_locations] = no_data_value

    del aster_b13_stdev_no_data_locations
    del aster_b14_stdev_no_data_locations

    # Create the ASTER EMIS standard deviation raster output tile
    logger.info('Creating an ASTER EMIS STDEV tile {}'.format(tile_name))
    util.Geo.generate_raster_file(gdal.GetDriverByName('GTiff'),
                                  tile_name,
                                  emis_stdev_data,
                                  samps, lines,
                                  transform,
                                  wkt,
                                  no_data_value,
                                  gdal.GDT_Float32)

    del emis_stdev_data


def generate_tiles(src_info, st_data_dir, url, wkt, no_data_value,
                   intermediate):
    """Generate tiles for emissivity standard deviation from ASTER data

    Args:
        src_info <SourceInfo>: Information about the source data
        st_data_dir <str>: Location of the ST data files
        url <str>: URL to retrieve the file from
        wkt <str>: Well-Known-Text describing the projection
        no_data_value <float>: Value to use for fill
        intermediate <bool>: Keep any intermediate products generated

    Returns:
        list(<str>): Standard deviation emissivity tile names
    """

    '''
    Process through the latitude and longitude ASTER tiles which cover
    the Landsat scene we are processing
    - Download them
    - Extract the Emissivity standard deviation bands 13 and 14
    - Generate the Landsat EMIS standard deviation from the 13 and 14 stdev
    - band data
    '''

    logger = logging.getLogger(__name__)

    # Read the ASTER GED tile list
    ged_tile_file = 'aster_ged_tile_list.txt'
    with open(os.path.join(st_data_dir, ged_tile_file)) as ged_file:
        tiles = [os.path.splitext(line.rstrip('\n'))[0] for line in ged_file]

    ls_emis_stdev_filenames = list()
    for (lat, lon) in [(lat, lon)
                       for lat in xrange(int(src_info.bound.south),
                                         int(src_info.bound.north)+1)
                       for lon in xrange(int(src_info.bound.west),
                                         int(src_info.bound.east)+1)]:

        # Build the base filename using the correct format
        filename = ''
        if lon < 0:
            filename = ASTER_GED_N_FORMAT.format(lat, lon)
        else:
            filename = ASTER_GED_P_FORMAT.format(lat, lon)

        # Skip the tile if it isn't in the ASTER GED database
        if filename not in tiles:
            logger.info('Skipping tile {} not in ASTER GED'.format(filename))
            continue

        # Build the output tile names
        ls_emis_stdev_tile_name = ''.join([filename, '_emis_stdev.tif'])

        # Read the ASTER data
        (aster_b13_stdev_data, aster_b14_stdev_data, samps, lines, transform,
         aster_data_available) = (
             extract_aster_data(url=url,
                                filename=filename,
                                intermediate=intermediate))

        # Fail if a tile can't be read, but it is in the ASTER GED
        if not aster_data_available:
            raise InaccessibleTileError(
                'Cannot reach tile {} in ASTER GED'.format(filename))

        # Add the tile names to the list for mosaic building and warping
        ls_emis_stdev_filenames.append(ls_emis_stdev_tile_name)

        generate_emis_stdev_tile(tile_name=ls_emis_stdev_tile_name,
                                 aster_b13_stdev_data=aster_b13_stdev_data,
                                 aster_b14_stdev_data=aster_b14_stdev_data,
                                 samps=samps,
                                 lines=lines,
                                 transform=transform,
                                 wkt=wkt,
                                 no_data_value=no_data_value)

        del aster_b13_stdev_data
        del aster_b14_stdev_data

    return ls_emis_stdev_filenames


def build_ls_emis_data(server_name, server_path, st_data_dir, src_info,
                       ls_emis_stdev_warped_name, no_data_value, intermediate):
    """Build estimated Landsat Emissivity Data

    Args:
        server_name <str>: Name of the ASTER GED server
        server_path <str>: Path on the ASTER GED server
        st_data_dir <str>: Location of the ST data files
        src_info <SourceInfo>: Information about the source data
        ls_emis_stdev_warped_name <str>: Path to warped emissivity stdev file
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
    # Save the source proj4 string to use during warping
    src_proj4 = ds_srs.ExportToProj4()

    (ls_emis_stdev_filenames) = (
        generate_tiles(src_info=src_info,
                       st_data_dir=st_data_dir,
                       url=url,
                       wkt=geographic_wkt,
                       no_data_value=no_data_value,
                       intermediate=intermediate))

    # Check to see that we downloaded at least one ASTER tile for processing.
    if len(ls_emis_stdev_filenames) == 0:
        raise NoTilesError('No ASTER tiles were downloaded')

    # Define the temporary names
    ls_emis_stdev_mosaic_name = 'landsat_emis_stdev_mosaic.tif'

    # Mosaic the estimated Landsat EMIS stdev tiles into the temp EMIS stdev
    logger.info('Building mosaic for estimated Landsat EMIS standard deviation')
    util.Geo.mosaic_tiles_into_one_raster(ls_emis_stdev_filenames,
                                          ls_emis_stdev_mosaic_name,
                                          no_data_value)

    if not intermediate:

        # Cleanup the estimated Landsat EMIS stdev tiles
        for emis_stdev_filename in ls_emis_stdev_filenames:
            if os.path.exists(emis_stdev_filename):
                os.unlink(emis_stdev_filename)

    # Warp estimated Landsat EMIS stdev to match the Landsat data
    logger.info('Warping estimated Landsat EMIS stdev to match Landsat data')
    emis_util.warp_raster(src_info, src_proj4, no_data_value,
                          ls_emis_stdev_mosaic_name, ls_emis_stdev_warped_name)

    if not intermediate:
        # Cleanup the temp files
        if os.path.exists(ls_emis_stdev_mosaic_name):
            os.unlink(ls_emis_stdev_mosaic_name)


def extract_warped_data(ls_emis_stdev_warped_name, no_data_value, intermediate):
    """Retrieves the warped image data

    Args:
        ls_emis_stdev_warped_name <str>: Path to warped emissivity stdev file
        no_data_value <float>: Value to use for fill
        intermediate <bool>: Keep any intermediate products generated

    Returns:
        <numpy.2darray>: Emissivity standard deviation data
        list(<int>): Emissivity stdev locations containing no data (fill) values
    """

    # Load the warped estimated Landsat EMIS stdev into memory
    ls_emis_stdev_data = emis_util.extract_raster_data(
        ls_emis_stdev_warped_name, 1)
    ls_emis_stdev_no_data_locations \
        = np.where(ls_emis_stdev_data == no_data_value)

    if not intermediate:
        # Cleanup the intermediate files since we have them in memory
        if os.path.exists(ls_emis_stdev_warped_name):
            os.unlink(ls_emis_stdev_warped_name)

    return (ls_emis_stdev_data, ls_emis_stdev_no_data_locations)


def generate_emissivity_data(xml_filename, server_name, server_path,
                             st_data_dir, no_data_value, intermediate):
    """Provides the main processing algorithm for generating the estimated
       Landsat emissivity product.  It produces the final emissivity product.

    Args:
        xml_filename <str>: Filename for the ESPA Metadata XML
        server_name <str>: Name of the ASTER GED server
        server_path <str>: Path on the ASTER GED server
        st_data_dir <str>: Location of the ST data files
        no_data_value <int>: No data (fill) value to use
        intermediate <bool>: Keep any intermediate products generated
    """

    logger = logging.getLogger(__name__)

    # XML metadata
    espa_metadata = Metadata(xml_filename)
    espa_metadata.parse()

    src_info = emis_util.retrieve_metadata_information(espa_metadata)

    # Determine output information
    sensor_code = emis_util.get_satellite_sensor_code(xml_filename)
    dataset = gdal.Open(src_info.toa.red.name)
    output_srs = osr.SpatialReference()
    output_srs.ImportFromWkt(dataset.GetProjection())
    output_transform = dataset.GetGeoTransform()
    samps = dataset.RasterXSize
    lines = dataset.RasterYSize
    del dataset

    ls_emis_stdev_warped_name = 'landsat_emis_stdev_warped.tif'

    # Build the estimated Landsat EMIS data from the ASTER GED data and
    # warp it to the Landsat scenes projection and image extents
    # For convenience the ASTER NDVI is also extracted and warped to the
    # Landsat scenes projection and image extents
    logger.info('Build thermal emissivity standard deviation band ')
    build_ls_emis_data(server_name=server_name,
                       server_path=server_path,
                       st_data_dir=st_data_dir,
                       src_info=src_info,
                       ls_emis_stdev_warped_name=ls_emis_stdev_warped_name,
                       no_data_value=no_data_value,
                       intermediate=intermediate)

    (ls_emis_stdev_data, ls_emis_stdev_no_data_locations) = (
        extract_warped_data(
            ls_emis_stdev_warped_name=ls_emis_stdev_warped_name,
            no_data_value=no_data_value,
            intermediate=intermediate))

    # Add the fill back into the results, since the may have been lost
    logger.info('Adding fill back into the estimated Landsat emissivity'
                ' stdev results')
    ls_emis_stdev_data[ls_emis_stdev_no_data_locations] = no_data_value

    # Memory cleanup
    del ls_emis_stdev_no_data_locations

    # Write emissivity standard deviation data and metadata
    ls_emis_stdev_img_filename = ''.join([xml_filename.split('.xml')[0],
                                          '_emis_stdev', '.img'])

    emis_util.write_emissivity_product(samps=samps,
                                       lines=lines,
                                       transform=output_transform,
                                       wkt=output_srs.ExportToWkt(),
                                       no_data_value=no_data_value,
                                       filename=ls_emis_stdev_img_filename,
                                       file_data=ls_emis_stdev_data)

    emis_util.add_emissivity_band_to_xml(espa_metadata=espa_metadata,
                                         filename=ls_emis_stdev_img_filename,
                                         sensor_code=sensor_code,
                                         no_data_value=no_data_value,
                                         band_type='stdev')

    # Memory cleanup
    del ls_emis_stdev_data


# Specify the no data value we will be using, it also matches the
# no_data_value for the ASTER data we extract and use
NO_DATA_VALUE = -9999


def main():
    """Generate Landsat EMIS and ASTER NDVI from ASTER GED tiles for the
       specified Landsat scene.
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

    logger.info('*** Begin Generate Estimated Landsat Emissivity ***')

    try:
        # Register all the gdal drivers
        gdal.AllRegister()

        # Get the data directory from the environment
        st_data_dir = emis_util.get_env_var('ST_DATA_DIR', None)

        # Call the main processing routine
        generate_emissivity_data(xml_filename=args.xml_filename,
                                 server_name=args.aster_ged_server_name,
                                 server_path=args.aster_ged_server_path,
                                 st_data_dir=st_data_dir,
                                 no_data_value=NO_DATA_VALUE,
                                 intermediate=args.intermediate)
    except Exception:
        logger.exception('Processing failed')
        sys.exit(1)  # EXIT FAILURE

    logger.info('*** Generate Estimated Landsat Emissivity - Complete ***')


if __name__ == '__main__':
    main()
