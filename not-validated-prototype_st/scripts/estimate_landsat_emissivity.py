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
import math
import datetime
from argparse import ArgumentParser
from collections import namedtuple


import requests
import numpy as np
from lxml import objectify as objectify
from osgeo import gdal, osr


from espa import Metadata
from st_exceptions import MissingBandError, NoTilesError


# Import local modules
import st_utilities as util


def extract_raster_data(name, band_number):
    """Extracts raster data for the specified dataset and band number

    Args:
        name <str>: Full path dataset name
        band_number <int>: Band number for the raster to extract

    Returns:
        <raster>: 2D raster array data
    """

    dataset = gdal.Open(name)
    if dataset is None:
        raise RuntimeError('GDAL failed to open {0}'.format(name))

    raster = (dataset.GetRasterBand(band_number)
              .ReadAsArray(0, 0, dataset.RasterXSize, dataset.RasterYSize))

    return raster


def data_resolution_and_size(name, x_min, x_max, y_min, y_max):
    """Calculates dataset resolution and retrieves size information

    Args:
        name <str>: Full path dataset name
        x_min <float>: Minimum longitude value
        x_max <float>: Maximum longitude value
        y_min <float>: Minimum latitude value
        y_max <float>: Maximum latitude value

    Returns:
        <float>: Longitude resolution
        <float>: Latitude resolution
        <int>: Samples in the data
        <int>: Lines in the data
    """

    dataset = gdal.Open(name)
    if dataset is None:
        raise RuntimeError('GDAL failed to open {0}'.format(name))

    return ((x_max - x_min) / float(dataset.RasterXSize),
            (y_max - y_min) / float(dataset.RasterYSize),
            dataset.RasterXSize, dataset.RasterYSize)


XYInfo = namedtuple('XYInfo',
                    ('x', 'y'))
BandInfo = namedtuple('BandInfo',
                      ('name', 'scale_factor', 'pixel_size'))
ToaInfo = namedtuple('ToaInfo',
                     ('green', 'red', 'nir', 'swir1', 'bt'))
ExtentInfo = namedtuple('ExtentInfo',
                        ('min', 'max'))
BoundInfo = namedtuple('BoundInfo',
                       ('north', 'south', 'east', 'west'))
SourceInfo = namedtuple('SourceInfo',
                        ('bound', 'extent', 'proj4', 'toa'))


def get_band_info(band):
    """Returns a populated BandInfo

    Args:
        band <xml_object>: Current band being processed

    Returns:
        <BandInfo>: Populated with band information
    """

    return BandInfo(name=str(band.file_name),
                    scale_factor=float(band.get('scale_factor')),
                    pixel_size=XYInfo(x=float(band.pixel_size.get('x')),
                                      y=float(band.pixel_size.get('y'))))


def extent_info(espa_metadata, band_info):
    """Returns a populated ExtentInfo

    Args:
        espa_metadata <espa.Metadata>: XML metadata

    Returns:
        <ExtentInfo>: Populated with extent information
    """

    # Determine the UTM projection corner points
    for corner_point in (espa_metadata.xml_object
                         .global_metadata
                         .projection_information.corner_point):
        if corner_point.get('location') == 'UL':
            extent_min_x = float(corner_point.get('x'))
            extent_max_y = float(corner_point.get('y'))
        if corner_point.get('location') == 'LR':
            extent_max_x = float(corner_point.get('x'))
            extent_min_y = float(corner_point.get('y'))

    '''
    Adjust the UTM coordinates for image extents becuse they are in
    center of pixel, and we need to supply the warping with actual
    extents
    '''
    extent_min_x = (extent_min_x - band_info.pixel_size.x * 0.5)
    extent_max_x = (extent_max_x + band_info.pixel_size.x * 0.5)
    extent_min_y = (extent_min_y - band_info.pixel_size.y * 0.5)
    extent_max_y = (extent_max_y + band_info.pixel_size.y * 0.5)

    return ExtentInfo(min=XYInfo(x=extent_min_x, y=extent_min_y),
                      max=XYInfo(x=extent_max_x, y=extent_max_y))


def bound_info(espa_metadata):
    """Returns a populated BoundInfo

    Args:
        espa_metadata <espa.Metadata>: XML metadata

    Returns:
        <BoundInfo>: Populated with boundary information
    """

    return BoundInfo(north=math.ceil(espa_metadata.xml_object.global_metadata
                                     .bounding_coordinates.north),
                     south=math.floor(espa_metadata.xml_object.global_metadata
                                      .bounding_coordinates.south),
                     east=math.ceil(espa_metadata.xml_object.global_metadata
                                    .bounding_coordinates.east),
                     west=math.floor(espa_metadata.xml_object.global_metadata
                                     .bounding_coordinates.west))


def retrieve_metadata_information(espa_metadata):
    """Reads required information from the metadata XML file

    Args:
        espa_metadata <espa.Metadata>: XML metadata

    Returns:
        <SourceInfo>: Populated with source information
    """

    bi_green = None
    bi_red = None
    bi_nir = None
    bi_swir1 = None
    bi_bt = None

    satellite = espa_metadata.xml_object.global_metadata.satellite

    # Find the TOA bands to extract information from
    for band in espa_metadata.xml_object.bands.band:
        if (band.get('product') == 'toa_refl' and
                band.get('name') == 'toa_band2'):
            bi_green = get_band_info(band)

        if (band.get('product') == 'toa_refl' and
                band.get('name') == 'toa_band3'):
            bi_red = get_band_info(band)

        if (band.get('product') == 'toa_refl' and
                band.get('name') == 'toa_band4'):
            bi_nir = get_band_info(band)

        if (band.get('product') == 'toa_refl' and
                band.get('name') == 'toa_band5'):
            bi_swir1 = get_band_info(band)

        if satellite == 'LANDSAT_8':
            if (band.get('product') == 'toa_bt' and
                    band.get('name') == 'bt_band11'):
                bi_bt = get_band_info(band)

                # Get the output proj4 string
                proj4 = util.Geo.get_proj4_projection_string(bi_bt.name)
        else:
            if (band.get('product') == 'toa_bt' and
                    band.get('name') == 'bt_band6'):
                bi_bt = get_band_info(band)

                # Get the output proj4 string
                proj4 = util.Geo.get_proj4_projection_string(bi_bt.name)

    # Error if we didn't find the required TOA bands in the data
    if bi_green is None:
        raise MissingBandError('Failed to find the TOA GREEN band'
                               ' in the input data')
    if bi_red is None:
        raise MissingBandError('Failed to find the TOA RED band'
                               ' in the input data')
    if bi_nir is None:
        raise MissingBandError('Failed to find the TOA NIR band'
                               ' in the input data')
    if bi_swir1 is None:
        raise MissingBandError('Failed to find the TOA SWIR1 band'
                               ' in the input data')
    if bi_bt is None:
        raise MissingBandError('Failed to find the TOA BT band'
                               ' in the input data')


    return SourceInfo(bound=bound_info(espa_metadata),
                      extent=extent_info(espa_metadata, bi_bt),
                      proj4=proj4,
                      toa=ToaInfo(green=bi_green,
                                  red=bi_red,
                                  nir=bi_nir,
                                  swir1=bi_swir1,
                                  bt=bi_bt))


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
    nir_data = extract_raster_data(src_info.toa.nir.name, 1)
    nir_no_data_locations = np.where(nir_data == no_data_value)
    nir_data = nir_data * src_info.toa.nir.scale_factor

    # RED ----------------------------------------------------------------
    red_data = extract_raster_data(src_info.toa.red.name, 1)
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
    green_data = extract_raster_data(src_info.toa.green.name, 1)
    green_no_data_locations = np.where(green_data == no_data_value)
    green_data = green_data * src_info.toa.green.scale_factor

    # SWIR1 --------------------------------------------------------------
    swir1_data = extract_raster_data(src_info.toa.swir1.name, 1)
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


def download_aster_ged_tile(url, h5_file_path):
    """Retrieves the specified tile from the host

    Args:
        url <str>: URL to retrieve the file from
        h5_file_path <str>: Full path on the remote system

    Raises:
        Exception: If issue transfering data
    """

    # Build the complete URL and download the tile
    url_path = ''.join([url, h5_file_path])
    status_code = util.Web.http_transfer_file(url_path, h5_file_path)

    # Check for and handle tiles that are not available in the
    # ASTER data
    if status_code != requests.codes['ok']:
        if status_code != requests.codes['not_found']:
            raise Exception('HTTP - Transfer Failed')


def extract_aster_data(url, filename, intermediate):
    """Extracts the internal band(s) data for later processing

    Args:
        url <str>: URL to retrieve the file from
        filename <str>: Base HDF filename to extract from
        intermediate <bool>: Keep any intermediate products generated

    Returns:
        <numpy.2darray>: Mean Band 13 data
        <numpy.2darray>: Mean Band 14 data
        <numpy.2darray>: SDev Band 13 data
        <numpy.2darray>: SDev Band 14 data
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

    download_aster_ged_tile(url=url, h5_file_path=h5_file_path)

    # There are cases where the emissivity data will not be available
    # (for example, in water regions).
    aster_b13_data = []
    aster_b14_data = []
    aster_b13_sdev_data = []
    aster_b14_sdev_data = []
    aster_ndvi_data = []
    samps = 0
    lines = 0
    geo_transform = []
    if not os.path.exists(h5_file_path):
        # The ASTER tile is not available, so don't try to process it
        return (aster_b13_data, aster_b14_data, aster_b13_sdev_data,
                aster_b14_sdev_data, aster_ndvi_data, samps, lines,
                geo_transform, False)

    # Define the sub-dataset names
    emis_ds_name = ''.join(['HDF5:"', h5_file_path,
                            '"://Emissivity/Mean'])
    emis_sdev_ds_name = ''.join(['HDF5:"', h5_file_path,
                            '"://Emissivity/SDev'])
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

    aster_b13_data = extract_raster_data(emis_ds_name, 4)
    aster_b14_data = extract_raster_data(emis_ds_name, 5)
    aster_b13_sdev_data = extract_raster_data(emis_sdev_ds_name, 4)
    aster_b14_sdev_data = extract_raster_data(emis_sdev_ds_name, 5)
    aster_ndvi_data = extract_raster_data(ndvi_ds_name, 1)
    aster_lat_data = extract_raster_data(lat_ds_name, 1)
    aster_lon_data = extract_raster_data(lon_ds_name, 1)

    # Determine the minimum and maximum latitude and longitude
    x_min = aster_lon_data.min()
    x_max = aster_lon_data.max()
    y_min = aster_lat_data.min()
    y_max = aster_lat_data.max()

    del aster_lon_data
    del aster_lat_data

    # Determine the resolution and dimensions of the ASTER data
    (x_res, y_res, samps, lines) = (
        data_resolution_and_size(lat_ds_name,
                                 x_min, x_max, y_min, y_max))

    # Remove the HDF5 tile since we no longer need it
    if not intermediate:
        if os.path.exists(h5_file_path):
            os.unlink(h5_file_path)

    # Build the geo transform
    geo_transform = [x_min, x_res, 0, y_max, 0, -y_res]

    return (aster_b13_data, aster_b14_data, aster_b13_sdev_data,
            aster_b14_sdev_data, aster_ndvi_data, samps, lines, geo_transform,
            True)


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

    # Save the no data and gap locations.
    aster_b13_stdev_gap_locations = np.where(aster_b13_stdev_data == 0)
    aster_b13_stdev_no_data_locations = np.where(aster_b13_stdev_data ==
        no_data_value)

    # Scale the data
    aster_b13_stdev_data = aster_b13_stdev_data * 0.0001

    # Save the no data and gap locations.
    aster_b14_stdev_gap_locations = np.where(aster_b14_stdev_data == 0)
    aster_b14_stdev_no_data_locations = np.where(aster_b14_stdev_data ==
        no_data_value)

    # Scale the data
    aster_b14_stdev_data = aster_b14_stdev_data * 0.0001

    # Create the estimated Landsat EMIS stdev data.
    emis_stdev_data = np.sqrt((aster_b13_stdev_data**2
        + aster_b14_stdev_data**2)/2)

    # Re-apply the no data and gap locations.
    emis_stdev_data[aster_b13_stdev_gap_locations] = 0
    emis_stdev_data[aster_b13_stdev_no_data_locations] = no_data_value
    emis_stdev_data[aster_b14_stdev_gap_locations] = 0
    emis_stdev_data[aster_b14_stdev_no_data_locations] = no_data_value

    del aster_b13_stdev_gap_locations
    del aster_b14_stdev_gap_locations
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
    """Generate tiles for emissivity mean, emissivity standard deviation,
       and NDVI from ASTER data

    Args:
        src_info <SourceInfo>: Information about the source data
        coefficients <CoefficientInfo>: coefficients for the math
        url <str>: URL to retrieve the file from
        wkt <str>: Well-Known-Text describing the projection
        no_data_value <float>: Value to use for fill
        intermediate <bool>: Keep any intermediate products generated

    Returns:
        list(<str>): Mean emissivity tile names
        list(<str>): Standard deviation emissivity tile names
        list(<str>): Mean ASTER NDVI tile names
    """

    '''
    Process through the latitude and longitude ASTER tiles which cover
    the Landsat scene we are processing
    - Download them
    - Extract the Emissivity mean bands 13 and 14
    - Extract the Emissivity standard deviation bands 13 and 14
    - Extract the NDVI
    - Generate the Landsat EMIS from the 13 and 14 band data
    - Generate the Landsat EMIS standard deviation from the 13 and 14 stdev
    - band data
    '''
    ls_emis_mean_filenames = list()
    ls_emis_stdev_filenames = list()
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
        ls_emis_stdev_tile_name = ''.join([filename, '_emis_stdev.tif'])
        aster_ndvi_tile_name = ''.join([filename, '_ndvi.tif'])

        # Read the ASTER data
        (aster_b13_data, aster_b14_data, aster_b13_stdev_data,
         aster_b14_stdev_data, aster_ndvi_data, samps, lines, transform,
         aster_data_available) = (
             extract_aster_data(url=url,
                                filename=filename,
                                intermediate=intermediate))

        # Skip the tile if it isn't available
        if not aster_data_available:
            continue

        # Add the tile names to the list for mosaic building and warping
        ls_emis_mean_filenames.append(ls_emis_tile_name)
        ls_emis_stdev_filenames.append(ls_emis_stdev_tile_name)
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

        generate_aster_ndvi_tile(tile_name=aster_ndvi_tile_name,
                                 ndvi_data=aster_ndvi_data,
                                 samps=samps,
                                 lines=lines,
                                 transform=transform,
                                 wkt=wkt,
                                 no_data_value=no_data_value)

        del aster_ndvi_data

    return (ls_emis_mean_filenames, ls_emis_stdev_filenames, 
            aster_ndvi_mean_filenames)


def warp_raster(target_info, src_proj4, no_data_value, src_name, dest_name):
    """Executes gdalwarp using the supplied information to warp to a specfic
       location and extent

    Args:
        src_info <float>: Destination pixel size for the x dimension
        src_proj4 <float>: Proj4 projection information for the source data
        src_name <str>: Name of the source data file
        dest_name <str>: Name of the destination data file
    """

    logger = logging.getLogger(__name__)

    cmd = ['gdalwarp', '-wm', '2048', '-multi',
           '-tr', str(target_info.toa.bt.pixel_size.x),
           str(target_info.toa.bt.pixel_size.y),
           '-s_srs', ''.join(["'", src_proj4, "'"]),
           '-t_srs', ''.join(["'", target_info.proj4, "'"]),
           '-of', 'GTiff',
           '-overwrite', '-te',
           str(target_info.extent.min.x), str(target_info.extent.min.y),
           str(target_info.extent.max.x), str(target_info.extent.max.y),
           '-srcnodata', str(no_data_value),
           '-dstnodata', str(no_data_value)]
    cmd.append(src_name)
    cmd.append(dest_name)

    cmd = ' '.join(cmd)

    output = ''
    try:
        logger.info('Executing [{0}]'.format(cmd))
        output = util.System.execute_cmd(cmd)
    except Exception:
        logger.error('Failed during warping')
        raise
    finally:
        if len(output) > 0:
            logger.info(output)


def build_ls_emis_data(server_name, server_path, src_info, coefficients,
                       ls_emis_warped_name, ls_emis_stdev_warped_name,
                       aster_ndvi_warped_name, no_data_value, intermediate):
    """Build estimated Landsat Emissivity Data

    Args:
        server_name <str>: Name of the ASTER GED server
        server_path <str>: Path on the ASTER GED server
        src_info <SourceInfo>: Information about the source data
        coefficients <CoefficientInfo>: coefficients for the math
        ls_emis_warped_name <str>: Path to the warped emissivity file
        ls_emis_stdev_warped_name <str>: Path to warped emissivity stdev file
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

    (ls_emis_mean_filenames, ls_emis_stdev_filenames, 
     aster_ndvi_mean_filenames) = (
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
    ls_emis_stdev_mosaic_name = 'landsat_emis_stdev_mosaic.tif'
    aster_ndvi_mosaic_name = 'aster_ndvi_mosaic.tif'

    # Mosaic the estimated Landsat EMIS tiles into the temp EMIS
    logger.info('Building mosaic for estimated Landsat EMIS')
    util.Geo.mosaic_tiles_into_one_raster(ls_emis_mean_filenames,
                                          ls_emis_mosaic_name,
                                          no_data_value)

    # Mosaic the estimated Landsat EMIS stdev tiles into the temp EMIS stdev
    logger.info('Building mosaic for estimated Landsat EMIS standard deviation')
    util.Geo.mosaic_tiles_into_one_raster(ls_emis_stdev_filenames,
                                          ls_emis_stdev_mosaic_name,
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

        # Cleanup the estimated Landsat EMIS stdev tiles
        for emis_stdev_filename in ls_emis_stdev_filenames:
            if os.path.exists(emis_stdev_filename):
                os.unlink(emis_stdev_filename)

        # Cleanup the ASTER NDVI tiles
        for ndvi_filename in aster_ndvi_mean_filenames:
            if os.path.exists(ndvi_filename):
                os.unlink(ndvi_filename)

    # Warp estimated Landsat EMIS to match the Landsat data
    logger.info('Warping estimated Landsat EMIS to match Landsat data')
    warp_raster(src_info, src_proj4, no_data_value,
                ls_emis_mosaic_name, ls_emis_warped_name)

    # Warp estimated Landsat EMIS stdev to match the Landsat data
    logger.info('Warping estimated Landsat EMIS stdev to match Landsat data')
    warp_raster(src_info, src_proj4, no_data_value,
                ls_emis_stdev_mosaic_name, ls_emis_stdev_warped_name)


    # Warp ASTER NDVI to match the Landsat data
    logger.info('Warping ASTER NDVI to match Landsat data')
    warp_raster(src_info, src_proj4, no_data_value,
                aster_ndvi_mosaic_name, aster_ndvi_warped_name)

    if not intermediate:
        # Cleanup the temp files
        if os.path.exists(ls_emis_mosaic_name):
            os.unlink(ls_emis_mosaic_name)
        if os.path.exists(ls_emis_stdev_mosaic_name):
            os.unlink(ls_emis_stdev_mosaic_name)
        if os.path.exists(aster_ndvi_mosaic_name):
            os.unlink(aster_ndvi_mosaic_name)


def extract_warped_data(ls_emis_warped_name, ls_emis_stdev_warped_name,
                        aster_ndvi_warped_name, no_data_value, intermediate):
    """Retrieves the warped image data with some massaging of ASTER NDVI

    Args:
        ls_emis_warped_name <str>: Path to the warped emissivity file
        ls_emis_stdev_warped_name <str>: Path to warped emissivity stdev file
        aster_ndvi_warped_name <str>: Path to the warped ASTER NDVI file
        no_data_value <float>: Value to use for fill
        intermediate <bool>: Keep any intermediate products generated

    Returns:
        <numpy.2darray>: Emissivity data
        list(<int>): Emissivity locations where gap data exists
        list(<int>): Emissivity locations containing no data (fill) values
        <numpy.2darray>: Emissivity standard deviation data
        list(<int>): Emissivity stdev deviation locations where gap data exists
        list(<int>): Emissivity stdev locations containing no data (fill) values
        <numpy.2darray>: ASTER NDVI data
        list(<int>): ASTER NDVI locations where gap data exists
        list(<int>): ASTER NDVI locations containing no data (fill) values
    """

    # Load the warped estimated Landsat EMIS into memory
    ls_emis_data = extract_raster_data(ls_emis_warped_name, 1)
    ls_emis_gap_locations = np.where(ls_emis_data == 0)
    ls_emis_no_data_locations = np.where(ls_emis_data == no_data_value)

    # Load the warped estimated Landsat EMIS stdev into memory
    ls_emis_stdev_data = extract_raster_data(ls_emis_stdev_warped_name, 1)
    ls_emis_stdev_gap_locations = np.where(ls_emis_stdev_data == 0)
    ls_emis_stdev_no_data_locations \
        = np.where(ls_emis_stdev_data == no_data_value)

    # Load the warped ASTER NDVI into memory
    aster_ndvi_data = extract_raster_data(aster_ndvi_warped_name, 1)
    aster_ndvi_gap_locations = np.where(aster_ndvi_data == 0)
    aster_ndvi_no_data_locations = np.where(aster_ndvi_data == no_data_value)

    # Turn all negative values to zero
    # Use a realy small value so that we don't have negative zero (-0.0)
    aster_ndvi_data[aster_ndvi_data < 0.0000001] = 0

    if not intermediate:
        # Cleanup the intermediate files since we have them in memory
        if os.path.exists(ls_emis_warped_name):
            os.unlink(ls_emis_warped_name)
        if os.path.exists(ls_emis_stdev_warped_name):
            os.unlink(ls_emis_stdev_warped_name)
        if os.path.exists(aster_ndvi_warped_name):
            os.unlink(aster_ndvi_warped_name)

    return (ls_emis_data, ls_emis_gap_locations, ls_emis_no_data_locations,
            ls_emis_stdev_data, ls_emis_stdev_gap_locations, 
            ls_emis_stdev_no_data_locations, aster_ndvi_data, 
            aster_ndvi_gap_locations, aster_ndvi_no_data_locations)


def write_emissivity_product(samps, lines, transform, wkt, no_data_value,
                             filename, file_data):
    """Creates the emissivity band file

    Args:
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
        filename <str>: Full path for the output file to create
        file_data <numpy.2darray>: Binary image data to be output
    """

    logger = logging.getLogger(__name__)

    logger.info('Creating {0}'.format(filename))
    util.Geo.generate_raster_file(gdal.GetDriverByName('ENVI'),
                                  filename,
                                  file_data,
                                  samps,
                                  lines,
                                  transform,
                                  wkt,
                                  no_data_value,
                                  gdal.GDT_Float32)

    hdr_filename = filename.replace('.img', '.hdr')
    logger.info('Updating {0}'.format(hdr_filename))
    util.Geo.update_envi_header(hdr_filename, no_data_value)

    # Remove the *.aux.xml file generated by GDAL
    aux_filename = filename.replace('.img', '.aux.xml')
    if os.path.exists(aux_filename):
        os.unlink(aux_filename)


def add_emissivity_band_to_xml(espa_metadata, filename, sensor_code,
                               no_data_value, band_type):
    """Adds the emissivity band to the Metadata XML file

    Args:
        espa_metadata <espa.Metadata>: XML metadata information
        filename <str>: Full path for the output file to create
        sensor_code <str>: Name prefix for the sensor 
        no_data_value <float>: Value to use for fill
        band_type <str>: Emissivity mean or standard deviation
    """

    logger = logging.getLogger(__name__)

    logger.info('Adding {0} to Metadata XML'.format(filename))

    # Create an element maker
    maker = objectify.ElementMaker(annotate=False, namespace=None, nsmap=None)

    source_product = 'toa_refl'

    # Find the TOA Band 1 to use for the specific band details
    base_band = None
    for band in espa_metadata.xml_object.bands.band:
        if (band.get('product') == source_product and
                band.get('name') == 'toa_band1'):
            base_band = band

    if base_band is None:
        raise MissingBandError('Failed to find the TOA BLUE band'
                               ' in the input data')

    emis_band = maker.band()
    emis_band.set('product', 'st_intermediate')
    emis_band.set('source', source_product)
    if band_type == 'mean':
        emis_band.set('name', 'emis')
    else: # band_type == 'stdev'
        emis_band.set('name', 'emis_stdev')
    emis_band.set('category', 'image')
    emis_band.set('data_type', 'FLOAT32')
    emis_band.set('nlines', base_band.attrib['nlines'])
    emis_band.set('nsamps', base_band.attrib['nsamps'])
    emis_band.set('fill_value', str(no_data_value))


    if band_type == 'mean':
        emis_band.short_name = maker.element('{0}EMIS'.format(sensor_code))
        emis_band.long_name = maker.element('Landsat emissivity estimated'
                                            ' from ASTER GED data')
    else: # band_type == 'stdev'
        emis_band.short_name \
            = maker.element('{0}EMIS_STDEV'.format(sensor_code))
        emis_band.long_name = maker.element('Landsat emissivity standard'
                                            ' deviation estimated from ASTER'
                                            ' GED data')
    emis_band.file_name = maker.element(filename)

    emis_band.pixel_size = base_band.pixel_size

    emis_band.resample_method = maker.element('none')
    emis_band.data_units = maker.element('Emissivity Coefficient')

    emis_band.valid_range = maker.element()
    emis_band.valid_range.set('min', '0.0')
    emis_band.valid_range.set('max', '1.0')

    emis_band.app_version = maker.element(util.Version.app_version())

    # Get the production date and time in string format
    # Strip the microseconds and add a Z
    date_now = ('{0}Z'.format(datetime.datetime.now()
                              .strftime('%Y-%m-%dT%H:%M:%S')))
    emis_band.production_date = maker.element(date_now)

    # Append the band to the XML
    espa_metadata.xml_object.bands.append(emis_band)

    # Validate the XML
    espa_metadata.validate()

    # Write it to the XML file
    espa_metadata.write()


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

    src_info = retrieve_metadata_information(espa_metadata)

    # Determine output information
    sensor_code = get_satellite_sensor_code(xml_filename)
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
    ls_emis_stdev_warped_name = 'landsat_emis_stdev_warped.tif'
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
                       ls_emis_stdev_warped_name=ls_emis_stdev_warped_name,
                       aster_ndvi_warped_name=aster_ndvi_warped_name,
                       no_data_value=no_data_value,
                       intermediate=intermediate)

    (ls_emis_data, ls_emis_gap_locations, ls_emis_no_data_locations,
     ls_emis_stdev_data, ls_emis_stdev_gap_locations,
     ls_emis_stdev_no_data_locations, aster_ndvi_data, aster_ndvi_gap_locations,
     aster_ndvi_no_data_locations) = (
         extract_warped_data(ls_emis_warped_name=ls_emis_warped_name,
                            ls_emis_stdev_warped_name=ls_emis_stdev_warped_name,
                            aster_ndvi_warped_name=aster_ndvi_warped_name,
                            no_data_value=no_data_value,
                            intermediate=intermediate))

    # Write emissivity standard deviation data and metadata
    ls_emis_stdev_img_filename = ''.join([xml_filename.split('.xml')[0],
                                    '_emis_stdev', '.img'])

    write_emissivity_product(samps=samps,
                             lines=lines,
                             transform=output_transform,
                             wkt=output_srs.ExportToWkt(),
                             no_data_value=no_data_value,
                             filename=ls_emis_stdev_img_filename,
                             file_data=ls_emis_stdev_data)

    add_emissivity_band_to_xml(espa_metadata=espa_metadata,
                               filename=ls_emis_stdev_img_filename,
                               sensor_code=sensor_code,
                               no_data_value=no_data_value,
                               band_type='stdev')

    # Memory cleanup
    del ls_emis_stdev_data

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

    write_emissivity_product(samps=samps,
                             lines=lines,
                             transform=output_transform,
                             wkt=output_srs.ExportToWkt(),
                             no_data_value=no_data_value,
                             filename=ls_emis_img_filename,
                             file_data=ls_emis_final)

    add_emissivity_band_to_xml(espa_metadata=espa_metadata,
                               filename=ls_emis_img_filename,
                               sensor_code=sensor_code,
                               no_data_value=no_data_value,
                               band_type='mean')

    # Memory cleanup
    del ls_emis_final


def retrieve_command_line_arguments():
    """Build the command line argument parser with some extra validation

    Returns:
        <args>: The command line arguments
    """

    description = ('Estimates Landsat Emissivity from ASTER GED data')
    parser = ArgumentParser(description=description)

    parser.add_argument('--version',
                        action='version',
                        version=util.Version.version_text())

    parser.add_argument('--xml',
                        action='store', dest='xml_filename',
                        required=False, default=None,
                        help='The XML metadata file to use')

    parser.add_argument('--aster-ged-server-name',
                        action='store', dest='aster_ged_server_name',
                        required=False, default=None,
                        help='Name of the ASTER GED server')

    parser.add_argument('--aster-ged-server-path',
                        action='store', dest='aster_ged_server_path',
                        required=False, default=None,
                        help='Path on the ASTER GED server')

    parser.add_argument('--intermediate',
                        action='store_true', dest='intermediate',
                        required=False, default=False,
                        help='Keep any intermediate products generated')

    parser.add_argument('--debug',
                        action='store_true', dest='debug',
                        required=False, default=False,
                        help='Turn debug messaging on')

    args = parser.parse_args()

    # Verify that the --xml parameter was specified
    if args.xml_filename is None:
        raise Exception('--xml must be specified on the command line')

    if args.aster_ged_server_name is None:
        raise Exception('--aster-ged-server-name must be specified on the'
                        ' command line')

    if args.aster_ged_server_name == '':
        raise Exception('The --aster-ged-server-name provided was empty')

    if args.aster_ged_server_path is None:
        raise Exception('--aster-ged-server-path must be specified on the'
                        ' command line')

    if args.aster_ged_server_path == '':
        raise Exception('The --aster-ged-server-path provided was empty')


    return args


def get_satellite_sensor_code(xml_filename):
    """Derives the satellite-sensor code from the XML filename

    Args:
        xml_filename <str>: Filename for the ESPA Metadata XML

    Returns:
        <str>: Satellite sensor code
    """

    old_prefixes = ['LT4', 'LT5', 'LE7', 'LT8', 'LC8', 'LO8']
    collection_prefixes = ['LT04', 'LT05', 'LE07', 'LT08', 'LC08', 'LO08']

    base_name = os.path.basename(xml_filename)

    satellite_sensor_code = base_name[0:3]
    if satellite_sensor_code in old_prefixes:
        return satellite_sensor_code

    satellite_sensor_code = base_name[0:4]
    if satellite_sensor_code in collection_prefixes:
        return satellite_sensor_code

    raise Exception('Satellite-Sensor code ({0}) not understood'
                    .format(satellite_sensor_code))


# Specify the no data value we will be using, it also matches the
# no_data_value for the ASTER data we extract and use
NO_DATA_VALUE = -9999


def main():
    """Generate Landsat EMIS and ASTER NDVI from ASTER GED tiles for the
       specified Landsat scene.
    """

    args = retrieve_command_line_arguments()

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
