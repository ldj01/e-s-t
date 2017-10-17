#! /usr/bin/env python

'''
    FILE: estimate_landsat_emissivity.py

    PURPOSE: Estimates a Landsat Emissivity product from ASTER Emissivity and
             NDVI.  The results are meant to be used for generation of a
             Surface Temperature product.

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
from collections import namedtuple


import numpy as np
from osgeo import gdal, osr


from espa import Metadata
from st_exceptions import NoTilesError


# Import local modules
import st_utilities as util
import emissivity_utilities as emis_util


CoefficientInfo = namedtuple('CoefficientInfo',
                             ('estimated_1', 'estimated_2', 'estimated_3',
                              'snow_emissivity', 'vegetation_coeff'))


def sensor_coefficients(satellite):
    """Determines the sensor specific coefficients

    Args:
        satellite <str>: Satellite we are currently processing

    Returns:
        <CoefficientInfo>: Populated coefficients
    """

    if satellite == 'LANDSAT_4':
        return CoefficientInfo(estimated_1=0.3222,
                               estimated_2=0.6498,
                               estimated_3=0.0272,
                               snow_emissivity=0.9883,
                               vegetation_coeff=0.9894)
    elif satellite == 'LANDSAT_5':
        return CoefficientInfo(estimated_1=-0.0723,
                               estimated_2=1.0521,
                               estimated_3=0.0195,
                               snow_emissivity=0.9856,
                               vegetation_coeff=0.9899)
    elif satellite == 'LANDSAT_7':
        return CoefficientInfo(estimated_1=0.2147,
                               estimated_2=0.7789,
                               estimated_3=0.0058,
                               snow_emissivity=0.9876,
                               vegetation_coeff=0.9896)
    elif satellite == 'LANDSAT_8':
        return CoefficientInfo(estimated_1=0.6820,
                               estimated_2=0.2578,
                               estimated_3=0.0584,
                               snow_emissivity=0.9904,
                               vegetation_coeff=0.9885)
    else:
        raise Exception('Unsupported satellite sensor')


def generate_landsat_ndvi(src_info, no_data_value):
    """Generate Landsat NDVI

    Args:
        src_info <SourceInfo>: Information about the source data
        no_data_value <int>: No data (fill) value to use

    Returns:
        <numpy.2darray>: Generated NDVI band data
        list(<int>): Locations containing no data (fill) values
    """

    logger = logging.getLogger(__name__)

    logger.info('Building TOA based NDVI band for Landsat data')

    # NIR ----------------------------------------------------------------
    nir_data = emis_util.extract_raster_data(src_info.toa.nir.name, 1)
    nir_no_data_locations = np.where(nir_data == no_data_value)
    nir_data = nir_data * src_info.toa.nir.scale_factor

    # RED ----------------------------------------------------------------
    red_data = emis_util.extract_raster_data(src_info.toa.red.name, 1)
    red_no_data_locations = np.where(red_data == no_data_value)
    red_data = red_data * src_info.toa.red.scale_factor

    # NDVI ---------------------------------------------------------------
    ndvi_data = ((nir_data - red_data) / (nir_data + red_data))

    # Cleanup no data locations
    ndvi_data[nir_no_data_locations] = no_data_value
    ndvi_data[red_no_data_locations] = no_data_value

    # Memory cleanup
    del red_data
    del nir_data
    del nir_no_data_locations
    del red_no_data_locations

    # Capture these before less than zero operation
    no_data_locations = np.where(ndvi_data == no_data_value)

    # Turn all negative values to zero
    # Use a realy small value so that we don't have negative zero (-0.0)
    ndvi_data[ndvi_data < 0.0000001] = 0

    return (ndvi_data, no_data_locations)


def snow_and_ndsi_locations(src_info, no_data_value):
    """Generate Landsat snow locations and NDSI nodata locations

    Args:
        src_info <SourceInfo>: Information about the source data
        no_data_value <int>: No data (fill) value to use

    Returns:
        list(<int>): Locations where we decided snow exists
        list(<int>): Locations containing no data (fill) values
    """

    logger = logging.getLogger(__name__)

    logger.info('Building TOA based NDSI band for Landsat data')

    # GREEN --------------------------------------------------------------
    green_data = emis_util.extract_raster_data(src_info.toa.green.name, 1)
    green_no_data_locations = np.where(green_data == no_data_value)
    green_data = green_data * src_info.toa.green.scale_factor

    # SWIR1 --------------------------------------------------------------
    swir1_data = emis_util.extract_raster_data(src_info.toa.swir1.name, 1)
    swir1_no_data_locations = np.where(swir1_data == no_data_value)
    swir1_data = swir1_data * src_info.toa.swir1.scale_factor

    # NDSI ---------------------------------------------------------------
    with np.errstate(divide='ignore'):
        ndsi_data = ((green_data - swir1_data) / (green_data + swir1_data))

    # Cleanup no data locations
    ndsi_data[green_no_data_locations] = no_data_value
    ndsi_data[swir1_no_data_locations] = no_data_value

    # Memory cleanup
    del green_data
    del swir1_data
    del green_no_data_locations
    del swir1_no_data_locations

    # Capture ndsi no data locations
    ndsi_no_data_locations = np.where(ndsi_data == no_data_value)

    # Save the locations for the specfied snow pixels
    logger.info('Determine snow pixel locations')
    snow_locations = np.where(ndsi_data > 0.4)

    # Memory cleanup
    del ndsi_data

    return (snow_locations, ndsi_no_data_locations)


ASTER_GED_N_FORMAT = 'AG100.v003.{0:02}.{1:04}.0001'
ASTER_GED_P_FORMAT = 'AG100.v003.{0:02}.{1:03}.0001'


def extract_aster_data(url, filename, intermediate):
    """Extracts the internal band(s) data for later processing

    Args:
        url <str>: URL to retrieve the file from
        filename <str>: Base HDF filename to extract from
        intermediate <bool>: Keep any intermediate products generated

    Returns:
        <numpy.2darray>: Mean Band 13 data
        <numpy.2darray>: Mean Band 14 data
        <numpy.2darray>: NDVI band data
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

    emis_util.download_aster_ged_tile(url=url, h5_file_path=h5_file_path)

    # There are cases where the emissivity data will not be available
    # (for example, in water regions).
    aster_b13_data = []
    aster_b14_data = []
    aster_ndvi_data = []
    samps = 0
    lines = 0
    geo_transform = []
    if not os.path.exists(h5_file_path):
        # The ASTER tile is not available, so don't try to process it
        return (aster_b13_data, aster_b14_data, aster_ndvi_data, samps, lines,
                geo_transform, False)

    # Define the sub-dataset names
    emis_ds_name = ''.join(['HDF5:"', h5_file_path,
                            '"://Emissivity/Mean'])
    ndvi_ds_name = ''.join(['HDF5:"', h5_file_path,
                            '"://NDVI/Mean'])
    lat_ds_name = ''.join(['HDF5:"', h5_file_path,
                           '"://Geolocation/Latitude'])
    lon_ds_name = ''.join(['HDF5:"', h5_file_path,
                           '"://Geolocation/Longitude'])

    logger.debug(emis_ds_name)
    logger.debug(ndvi_ds_name)
    logger.debug(lat_ds_name)
    logger.debug(lon_ds_name)

    aster_b13_data = emis_util.extract_raster_data(emis_ds_name, 4)
    aster_b14_data = emis_util.extract_raster_data(emis_ds_name, 5)
    aster_ndvi_data = emis_util.extract_raster_data(ndvi_ds_name, 1)
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

    return (aster_b13_data, aster_b14_data, aster_ndvi_data, samps, lines, 
            geo_transform, True)


def generate_estimated_emis_tile(coefficients, tile_name,
                                 aster_b13_data, aster_b14_data,
                                 samps, lines, transform,
                                 wkt, no_data_value):
    """Generate emissivity values for the tile

    Args:
        coefficients <CoefficientInfo>: coefficients for the math
        tile_name <str>: Filename to create for the tile
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

    # Save the no data and gap locations.
    aster_b13_gap_locations = np.where(aster_b13_data == 0)
    aster_b13_no_data_locations = np.where(aster_b13_data == no_data_value)

    # Scale the data
    aster_b13_data = aster_b13_data * 0.001

    # Save the no data and gap locations.
    aster_b14_gap_locations = np.where(aster_b14_data == 0)
    aster_b14_no_data_locations = np.where(aster_b14_data == no_data_value)

    # Scale the data
    aster_b14_data = aster_b14_data * 0.001

    # ------------------------------------------------------------
    # Create the estimated Landsat EMIS data
    emis_data = (coefficients.estimated_1 * aster_b13_data +
                 coefficients.estimated_2 * aster_b14_data +
                 coefficients.estimated_3)

    # Re-apply the no data and gap locations.
    emis_data[aster_b13_gap_locations] = 0
    emis_data[aster_b13_no_data_locations] = no_data_value
    emis_data[aster_b14_gap_locations] = 0
    emis_data[aster_b14_no_data_locations] = no_data_value

    del aster_b13_gap_locations
    del aster_b14_gap_locations
    del aster_b13_no_data_locations
    del aster_b14_no_data_locations

    # Create the estimated Landsat EMIS raster output tile
    logger.info('Creating an estimated Landsat EMIS tile {}'.format(tile_name))
    util.Geo.generate_raster_file(gdal.GetDriverByName('GTiff'),
                                  tile_name,
                                  emis_data,
                                  samps, lines,
                                  transform,
                                  wkt,
                                  no_data_value,
                                  gdal.GDT_Float32)

    del emis_data


def generate_aster_ndvi_tile(tile_name, ndvi_data,
                             samps, lines, transform,
                             wkt, no_data_value):
    """Generate ASTER NDVI values for the tile

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

    # Save the no data locations.
    ndvi_no_data_locations = np.where(ndvi_data == no_data_value)

    # Scale the data
    data = ndvi_data * 0.01

    # Re-apply the no data locations.
    data[ndvi_no_data_locations] = no_data_value

    # Create the ASTER NDVI raster output tile
    logger.info('Creating an ASTER NDVI tile {}'.format(tile_name))
    util.Geo.generate_raster_file(gdal.GetDriverByName('GTiff'),
                                  tile_name,
                                  data,
                                  samps, lines,
                                  transform,
                                  wkt,
                                  no_data_value,
                                  gdal.GDT_Float32)

    del data


def generate_tiles(src_info, coefficients, url, wkt,
                   no_data_value, intermediate):
    """Generate tiles for emissivity mean and NDVI from ASTER data

    Args:
        src_info <SourceInfo>: Information about the source data
        coefficients <CoefficientInfo>: coefficients for the math
        url <str>: URL to retrieve the file from
        wkt <str>: Well-Known-Text describing the projection
        no_data_value <float>: Value to use for fill
        intermediate <bool>: Keep any intermediate products generated

    Returns:
        list(<str>): Mean emissivity tile names
        list(<str>): Mean ASTER NDVI tile names
    """

    '''
    Process through the latitude and longitude ASTER tiles which cover
    the Landsat scene we are processing
    - Download them
    - Extract the Emissivity mean bands 13 and 14
    - Extract the NDVI
    - Generate the Landsat EMIS from the 13 and 14 band data
    '''
    ls_emis_mean_filenames = list()
    aster_ndvi_mean_filenames = list()
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

        # Build the output tile names
        ls_emis_tile_name = ''.join([filename, '_emis.tif'])
        aster_ndvi_tile_name = ''.join([filename, '_ndvi.tif'])

        # Read the ASTER data
        (aster_b13_data, aster_b14_data, aster_ndvi_data, samps, lines, 
         transform, aster_data_available) = (
             extract_aster_data(url=url,
                                filename=filename,
                                intermediate=intermediate))

        # Skip the tile if it isn't available
        if not aster_data_available:
            continue

        # Add the tile names to the list for mosaic building and warping
        ls_emis_mean_filenames.append(ls_emis_tile_name)
        aster_ndvi_mean_filenames.append(aster_ndvi_tile_name)

        generate_estimated_emis_tile(coefficients=coefficients,
                                     tile_name=ls_emis_tile_name,
                                     aster_b13_data=aster_b13_data,
                                     aster_b14_data=aster_b14_data,
                                     samps=samps,
                                     lines=lines,
                                     transform=transform,
                                     wkt=wkt,
                                     no_data_value=no_data_value)

        del aster_b13_data
        del aster_b14_data

        generate_aster_ndvi_tile(tile_name=aster_ndvi_tile_name,
                                 ndvi_data=aster_ndvi_data,
                                 samps=samps,
                                 lines=lines,
                                 transform=transform,
                                 wkt=wkt,
                                 no_data_value=no_data_value)

        del aster_ndvi_data

    return (ls_emis_mean_filenames, aster_ndvi_mean_filenames)


def build_ls_emis_data(server_name, server_path, src_info, coefficients,
                       ls_emis_warped_name, aster_ndvi_warped_name, 
                       no_data_value, intermediate):
    """Build estimated Landsat Emissivity Data

    Args:
        server_name <str>: Name of the ASTER GED server
        server_path <str>: Path on the ASTER GED server
        src_info <SourceInfo>: Information about the source data
        coefficients <CoefficientInfo>: coefficients for the math
        ls_emis_warped_name <str>: Path to the warped emissivity file
        aster_ndvi_warped_name <str>: Path to the warped ASTER NDVI file
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

    (ls_emis_mean_filenames, aster_ndvi_mean_filenames) = (
        generate_tiles(src_info=src_info,
                       coefficients=coefficients,
                       url=url,
                       wkt=geographic_wkt,
                       no_data_value=no_data_value,
                       intermediate=intermediate))

    # Check to see that we downloaded at least one ASTER tile for processing.
    if len(ls_emis_mean_filenames) == 0:
        raise NoTilesError('No ASTER tiles were downloaded')

    # Define the temporary names
    ls_emis_mosaic_name = 'landsat_emis_mosaic.tif'
    aster_ndvi_mosaic_name = 'aster_ndvi_mosaic.tif'

    # Mosaic the estimated Landsat EMIS tiles into the temp EMIS
    logger.info('Building mosaic for estimated Landsat EMIS')
    util.Geo.mosaic_tiles_into_one_raster(ls_emis_mean_filenames,
                                          ls_emis_mosaic_name,
                                          no_data_value)

    # Mosaic the ASTER NDVI tiles into the temp NDVI
    logger.info('Building mosaic for ASTER NDVI')
    util.Geo.mosaic_tiles_into_one_raster(aster_ndvi_mean_filenames,
                                          aster_ndvi_mosaic_name,
                                          no_data_value)

    if not intermediate:
        # Cleanup the estimated Landsat EMIS tiles
        for emis_filename in ls_emis_mean_filenames:
            if os.path.exists(emis_filename):
                os.unlink(emis_filename)

        # Cleanup the ASTER NDVI tiles
        for ndvi_filename in aster_ndvi_mean_filenames:
            if os.path.exists(ndvi_filename):
                os.unlink(ndvi_filename)

    # Warp estimated Landsat EMIS to match the Landsat data
    logger.info('Warping estimated Landsat EMIS to match Landsat data')
    emis_util.warp_raster(src_info, src_proj4, no_data_value,
                ls_emis_mosaic_name, ls_emis_warped_name)

    # Warp ASTER NDVI to match the Landsat data
    logger.info('Warping ASTER NDVI to match Landsat data')
    emis_util.warp_raster(src_info, src_proj4, no_data_value,
                aster_ndvi_mosaic_name, aster_ndvi_warped_name)

    if not intermediate:
        # Cleanup the temp files
        if os.path.exists(ls_emis_mosaic_name):
            os.unlink(ls_emis_mosaic_name)
        if os.path.exists(aster_ndvi_mosaic_name):
            os.unlink(aster_ndvi_mosaic_name)


def extract_warped_data(ls_emis_warped_name, aster_ndvi_warped_name, 
                        no_data_value, intermediate):
    """Retrieves the warped image data with some massaging of ASTER NDVI

    Args:
        ls_emis_warped_name <str>: Path to the warped emissivity file
        aster_ndvi_warped_name <str>: Path to the warped ASTER NDVI file
        no_data_value <float>: Value to use for fill
        intermediate <bool>: Keep any intermediate products generated

    Returns:
        <numpy.2darray>: Emissivity data
        list(<int>): Emissivity locations where gap data exists
        list(<int>): Emissivity locations containing no data (fill) values
        <numpy.2darray>: Emissivity standard deviation data
        <numpy.2darray>: ASTER NDVI data
        list(<int>): ASTER NDVI locations where gap data exists
        list(<int>): ASTER NDVI locations containing no data (fill) values
    """

    # Load the warped estimated Landsat EMIS into memory
    ls_emis_data = emis_util.extract_raster_data(ls_emis_warped_name, 1)
    ls_emis_gap_locations = np.where(ls_emis_data == 0)
    ls_emis_no_data_locations = np.where(ls_emis_data == no_data_value)

    # Load the warped ASTER NDVI into memory
    aster_ndvi_data = emis_util.extract_raster_data(aster_ndvi_warped_name, 1)
    aster_ndvi_gap_locations = np.where(aster_ndvi_data == 0)
    aster_ndvi_no_data_locations = np.where(aster_ndvi_data == no_data_value)

    # Turn all negative values to zero
    # Use a realy small value so that we don't have negative zero (-0.0)
    aster_ndvi_data[aster_ndvi_data < 0.0000001] = 0

    if not intermediate:
        # Cleanup the intermediate files since we have them in memory
        if os.path.exists(ls_emis_warped_name):
            os.unlink(ls_emis_warped_name)
        if os.path.exists(aster_ndvi_warped_name):
            os.unlink(aster_ndvi_warped_name)

    return (ls_emis_data, ls_emis_gap_locations, ls_emis_no_data_locations,
            aster_ndvi_data, aster_ndvi_gap_locations, 
            aster_ndvi_no_data_locations)


def generate_emissivity_data(xml_filename, server_name, server_path,
                             no_data_value, intermediate):
    """Provides the main processing algorithm for generating the estimated
       Landsat emissivity product.  It produces the final emissivity product.

    Args:
        xml_filename <str>: Filename for the ESPA Metadata XML
        server_name <str>: Name of the ASTER GED server
        server_path <str>: Path on the ASTER GED server
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

    # Initialize coefficients.
    ASTER_GED_WATER = 0.988
    coefficients = sensor_coefficients(espa_metadata.xml_object
                                       .global_metadata.satellite)

    # ====================================================================
    # Build NDVI in memory
    (ls_ndvi_data, ndvi_no_data_locations) = (
        generate_landsat_ndvi(src_info, no_data_value))

    if intermediate:
        logger.info('Writing Landsat NDVI raster')
        util.Geo.generate_raster_file(gdal.GetDriverByName('GTiff'),
                                      'internal_landsat_ndvi.tif',
                                      ls_ndvi_data,
                                      samps,
                                      lines,
                                      output_transform,
                                      output_srs.ExportToWkt(),
                                      no_data_value,
                                      gdal.GDT_Float32)

    # ====================================================================
    # Determine NDSI and Snow locations
    (snow_locations, ndsi_no_data_locations) = (
        snow_and_ndsi_locations(src_info, no_data_value))

    ls_emis_warped_name = 'landsat_emis_warped.tif'
    aster_ndvi_warped_name = 'aster_ndvi_warped.tif'

    # Build the estimated Landsat EMIS data from the ASTER GED data and
    # warp it to the Landsat scenes projection and image extents
    # For convenience the ASTER NDVI is also extracted and warped to the
    # Landsat scenes projection and image extents
    logger.info('Build thermal emissivity band and retrieve ASTER NDVI')
    build_ls_emis_data(server_name=server_name,
                       server_path=server_path,
                       src_info=src_info,
                       coefficients=coefficients,
                       ls_emis_warped_name=ls_emis_warped_name,
                       aster_ndvi_warped_name=aster_ndvi_warped_name,
                       no_data_value=no_data_value,
                       intermediate=intermediate)

    (ls_emis_data, ls_emis_gap_locations, ls_emis_no_data_locations,
     aster_ndvi_data, aster_ndvi_gap_locations, aster_ndvi_no_data_locations) \
         = (extract_warped_data(ls_emis_warped_name=ls_emis_warped_name,
                            aster_ndvi_warped_name=aster_ndvi_warped_name,
                            no_data_value=no_data_value,
                            intermediate=intermediate))

    # Find water locations using the value of water in ASTER GED
    water_locations = np.where(ls_emis_data > ASTER_GED_WATER)

    # Replace NDVI values greater than 1 with 1
    ls_ndvi_data[ls_ndvi_data > 1.0] = 1
    aster_ndvi_data[aster_ndvi_data > 1.0] = 1

    logger.info('Normalizing Landsat and ASTER NDVI')
    # Normalize Landsat NDVI by max value
    max_ls_ndvi = ls_ndvi_data.max()
    min_ls_ndvi = ls_ndvi_data.min()
    logger.info('Max LS NDVI {0}'.format(max_ls_ndvi))
    ls_ndvi_data = ls_ndvi_data / float(max_ls_ndvi)

    if intermediate:
        logger.info('Writing Landsat NDVI NORM MAX raster')
        util.Geo.generate_raster_file(gdal.GetDriverByName('GTiff'),
                                      'internal_landsat_ndvi_norm_max.tif',
                                      ls_ndvi_data,
                                      samps,
                                      lines,
                                      output_transform,
                                      output_srs.ExportToWkt(),
                                      no_data_value,
                                      gdal.GDT_Float32)

    # Normalize ASTER NDVI by max value
    max_aster_ndvi = aster_ndvi_data.max()
    logger.info('Max ASTER NDVI {0}'.format(max_aster_ndvi))
    aster_ndvi_data = aster_ndvi_data / float(max_aster_ndvi)

    if intermediate:
        logger.info('Writing Aster NDVI NORM MAX raster')
        util.Geo.generate_raster_file(gdal.GetDriverByName('GTiff'),
                                      'internal_aster_ndvi_norm_max.tif',
                                      aster_ndvi_data,
                                      samps,
                                      lines,
                                      output_transform,
                                      output_srs.ExportToWkt(),
                                      no_data_value,
                                      gdal.GDT_Float32)

    # Soil - From prototype code variable name
    logger.info('Calculating bare soil component')

    # Get pixels with significant bare soil component 
    bare_locations = np.where(aster_ndvi_data < 0.5)

    # Only calculate soil component for these pixels
    ls_emis_bare = ((ls_emis_data[bare_locations]
                     - 0.975 * aster_ndvi_data[bare_locations])
                    / (1 - aster_ndvi_data[bare_locations]))

    # Memory cleanup
    del aster_ndvi_data

    # Calculate veg adjustment with Landsat
    logger.info('Calculating EMIS Final')

    # Adjust estimated Landsat EMIS for vegetation and snow, to generate
    # the final Landsat EMIS data
    logger.info('Adjusting estimated EMIS for vegetation')
    ls_emis_final = (coefficients.vegetation_coeff * ls_ndvi_data +
                     ls_emis_data * (1.0 - ls_ndvi_data))

    # Calculate fractional vegetation cover
    fv_L = 1.0 - (max_ls_ndvi - ls_ndvi_data) / (max_ls_ndvi - min_ls_ndvi)

    # Memory cleanup
    del ls_emis_data
    del ls_ndvi_data

    # Add soil component pixels
    ls_emis_final[bare_locations] = ls_emis_bare

    # Memory cleanup
    del ls_emis_bare
    del bare_locations

    # Set fill values on granule edge to nan
    fill_locations = np.where(np.isnan(fv_L))
    ls_emis_final[fill_locations] = np.nan

    # Memory cleanup
    del fv_L
    del fill_locations

    # Final check for emissivity values greater than 1.  Reset values greater
    # than 1 to nominal veg/water value (should be very few, if any)
    ls_emis_final[ls_emis_final > 1.0] = ASTER_GED_WATER

    # Medium snow
    logger.info('Adjusting estimated EMIS for snow')
    ls_emis_final[snow_locations] = coefficients.snow_emissivity

    # Memory cleanup
    del snow_locations

    # Reset water values
    ls_emis_final[water_locations] = ASTER_GED_WATER 

    # Memory cleanup
    del water_locations

    # Add the fill and scan gaps and ASTER gaps back into the results,
    # since they may have been lost
    logger.info('Adding fill and data gaps back into the estimated'
                ' Landsat emissivity results')
    ls_emis_final[ls_emis_no_data_locations] = no_data_value
    ls_emis_final[ls_emis_gap_locations] = no_data_value
    ls_emis_final[aster_ndvi_no_data_locations] = no_data_value
    ls_emis_final[aster_ndvi_gap_locations] = no_data_value
    ls_emis_final[ndvi_no_data_locations] = no_data_value
    ls_emis_final[ndsi_no_data_locations] = no_data_value

    # Memory cleanup
    del ls_emis_no_data_locations
    del ls_emis_gap_locations
    del aster_ndvi_no_data_locations
    del aster_ndvi_gap_locations
    del ndvi_no_data_locations
    del ndsi_no_data_locations 

    # Write emissivity data and metadata
    ls_emis_img_filename = ''.join([xml_filename.split('.xml')[0],
                                    '_emis', '.img'])

    emis_util.write_emissivity_product(samps=samps,
                             lines=lines,
                             transform=output_transform,
                             wkt=output_srs.ExportToWkt(),
                             no_data_value=no_data_value,
                             filename=ls_emis_img_filename,
                             file_data=ls_emis_final)

    emis_util.add_emissivity_band_to_xml(espa_metadata=espa_metadata,
                               filename=ls_emis_img_filename,
                               sensor_code=sensor_code,
                               no_data_value=no_data_value,
                               band_type='mean')

    # Memory cleanup
    del ls_emis_final



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

        # Call the main processing routine
        generate_emissivity_data(xml_filename=args.xml_filename,
                                 server_name=args.aster_ged_server_name,
                                 server_path=args.aster_ged_server_path,
                                 no_data_value=NO_DATA_VALUE,
                                 intermediate=args.intermediate)
    except Exception:
        logger.exception('Processing failed')
        sys.exit(1)  # EXIT FAILURE

    logger.info('*** Generate Estimated Landsat Emissivity - Complete ***')


if __name__ == '__main__':
    main()
