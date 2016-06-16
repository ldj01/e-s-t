#! /usr/bin/env python

'''
    FILE: estimate_landsat_emissivity.py

    PURPOSE: Estimates a Landsat Emissivity product from ASTER Emissivity and
             NDVI.  The results are meant to be used for generation of a Land
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
from lst_exceptions import MissingBandError


# Import local modules
import lst_utilities as util


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
    """

    dataset = gdal.Open(name)
    if dataset is None:
        raise RuntimeError('GDAL failed to open {0}'.format(name))

    return ((x_max - x_min) / float(dataset.RasterXSize),
            (y_max - y_min) / float(dataset.RasterYSize),
            dataset.RasterXSize, dataset.RasterYSize)


def warp_raster(target_info, src_proj4, no_data_value, src_name, dest_name):
    """Executes gdalwarp using the supplied  information to warp to a specfic
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
        band <xml_object>: Current being processed

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

    return ExtentInfo(min=XYInfo(x=extent_min_x,
                                 y=extent_min_y),
                      max=XYInfo(x=extent_max_x,
                                 y=extent_max_y))


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

        if (band.get('product') == 'toa_bt' and
                band.get('category') == 'image'):
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
                             ('estimated_1', 'estimated_2', 'estimated_3', 'snow_emissivity',
                              'vegetation_coeff'))


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
        <data>: Generated NDVI band data
        <array>: Array of locations containing no data (fill) values
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
    """Generate Landsat NDVI

    Args:
        src_info <SourceInfo>: Information about the source data
        no_data_value <int>: No data (fill) value to use

    Returns:
        <array>: Array of locations where we decided snow exists
        <array>: Array of locations containing no data (fill) values
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


def generate_emissivity_data(espa_metadata):
    """Generate the Landsat Emissivity band data
    Args:
    Returns:
    """

    logger = logging.getLogger(__name__)


class EstimateLandsatEmissivity(object):
    '''
    Description:
        Defines the processor for generating the estimated Landsat Emissivity
        product.
    '''

    def __init__(self, xml_filename,
                 aster_ged_server_name, aster_ged_server_path,
                 keep_intermediate_data=False):
        super(EstimateLandsatEmissivity, self).__init__()

        # Keep local copies of these
        self.xml_filename = xml_filename
        self.keep_intermediate_data = keep_intermediate_data

        server_name = aster_ged_server_name
        server_path = aster_ged_server_path

        # Specify the no data value we will be using, it also matches the
        # no_data_value for the ASTER data we extract and use
        self.no_data_value = -9999

        # Specify the base URL to use for retrieving the ASTER GED data
        self.base_url_path = ''.join(['http://', server_name, server_path])

        # ASTER GED filename formats
        # AG100.v003.44.-077.0001.h5  (Negative Longitude)
        # or
        # AG100.v003.44.077.0001.h5   (Positive Longitude)
        # Setup the formats to apply to the URL for retrieving ASTER GED tiles
        self.file_n_format = 'AG100.v003.{0:02}.{1:04}.0001'
        self.file_p_format = 'AG100.v003.{0:02}.{1:03}.0001'

    def build_ls_emis_data(self, src_info, coefficients):
        '''
        Description:
            Download the ASTER GED tiles that encompass our Landsat scene and
            extract the bands required to generate the Landsat Emissivity data.
            Mosaic the Landsat Emissivity tiles together and then warp them to
            the projection and image extents of the Landsat scenes.  For
            convenience the ASTER NDVI is also extracted and warped to the
            Landsat scenes projection and image extents

        Returns:
            ls_emis_warped_name
              - The name of the reprojected Landsat Emissivity data.
            aster_ndvi_warped_name
              - The name of the reprojected ASTER NDVI data.
        '''

        logger = logging.getLogger(__name__)

        # The ASTER data is in geographic projection so specify that here
        ds_srs = osr.SpatialReference()
        ds_srs.ImportFromEPSG(4326)

        # Process through the lattitude and longitude ASTER tiles which cover
        # the Landsat scene we are processing
        # - Download them
        # - Extract the Emissivity bands 13 and 14 as well as the NDVI
        # - Generate the Landsat EMIS from the 13 and 14 band data
        ls_emis_mean_filenames = list()
        aster_ndvi_mean_filenames = list()
        for lat in xrange(int(src_info.bound.south), int(src_info.bound.north)+1):
            for lon in xrange(int(src_info.bound.west), int(src_info.bound.east)+1):
                # Build the base filename using the correct format
                filename = ''
                if lon < 0:
                    filename = self.file_n_format.format(lat, lon)
                else:
                    filename = self.file_p_format.format(lat, lon)

                # Build the HDF5 filename for the tile
                h5_file_path = ''.join([filename, '.h5'])

                # Build the URL and download the tile
                url_path = ''.join([self.base_url_path, h5_file_path])
                status_code = util.Web.http_transfer_file(url_path,
                                                          h5_file_path)

                # Check for and handle tiles that are not available in the
                # ASTER data
                if status_code != requests.codes['ok']:
                    if status_code != requests.codes['not_found']:
                        raise Exception('HTTP - Transfer Failed')
                    else:
                        # Advance to the next tile
                        continue

                # ------------------------------------------------------------
                # Build the output tile names
                ls_emis_tile_name = ''.join([filename, '_emis.tif'])
                aster_ndvi_tile_name = ''.join([filename, '_ndvi.tif'])

                # Add the tile names to the list for mosaic building and
                # warping
                ls_emis_mean_filenames.append(ls_emis_tile_name)
                aster_ndvi_mean_filenames.append(aster_ndvi_tile_name)

                # Get the sub-dataset name for Emissivity-Mean
                emis_ds_name = ''.join(['HDF5:"', h5_file_path,
                                        '"://Emissivity/Mean'])
                # Get the sub-dataset name for NDVI-ean
                ndvi_ds_name = ''.join(['HDF5:"', h5_file_path,
                                        '"://NDVI/Mean'])
                # Get the sub-dataset name for Lattitude
                lat_ds_name = ''.join(['HDF5:"', h5_file_path,
                                       '"://Geolocation/Latitude'])
                # Get the sub-dataset name for Longitude
                lon_ds_name = ''.join(['HDF5:"', h5_file_path,
                                       '"://Geolocation/Longitude'])

                logger.debug(emis_ds_name)
                logger.debug(ndvi_ds_name)
                logger.debug(lat_ds_name)
                logger.debug(lon_ds_name)

                # ------------------------------------------------------------
                try:
                    aster_b13_data = extract_raster_data(emis_ds_name, 4)
                    aster_b14_data = extract_raster_data(emis_ds_name, 5)
                    aster_ndvi_data = extract_raster_data(ndvi_ds_name, 1)
                    aster_lat_data = extract_raster_data(lat_ds_name, 1)
                    aster_lon_data = extract_raster_data(lon_ds_name, 1)

                except Exception:
                    logger.exception('Extracting ASTER data from tile')
                    raise

                # ------------------------------------------------------------
                # Determine the minimum and maximum latitude and longitude
                x_min = aster_lon_data.min()
                x_max = aster_lon_data.max()
                y_min = aster_lat_data.min()
                y_max = aster_lat_data.max()

                # Determine the resolution of the ASTER data
                (x_res, y_res, x_dim, y_dim) = (
                    data_resolution_and_size(lat_ds_name,
                                             x_min, x_max, y_min, y_max))

                # Remove the HDF5 tile since we have extracted all the info we
                # need from it
                if not self.keep_intermediate_data:
                    if os.path.exists(h5_file_path):
                        os.unlink(h5_file_path)

                # Build the geo transform to apply to the raster tile
                geo_transform = [x_min, x_res, 0, y_max, 0, -y_res]

                # ------------------------------------------------------------
                # Save the no data and gap locations.
                aster_b13_gap_locations = np.where(aster_b13_data == 0)
                aster_b13_no_data_locations = (
                    np.where(aster_b13_data == self.no_data_value))

                # Scale the data
                aster_b13_data = aster_b13_data * 0.001

                # ------------------------------------------------------------
                # Save the no data and gap locations.
                aster_b14_gap_locations = np.where(aster_b14_data == 0)
                aster_b14_no_data_locations = (
                    np.where(aster_b14_data == self.no_data_value))

                # Scale the data
                aster_b14_data = aster_b14_data * 0.001

                # ------------------------------------------------------------
                # Create the estimated Landsat EMIS data
                ls_emis_data = (coefficients.estimated_1 * aster_b13_data +
                                coefficients.estimated_2 * aster_b14_data +
                                coefficients.estimated_3)

                # Re-apply the no data and gap locations.
                ls_emis_data[aster_b13_gap_locations] = 0
                ls_emis_data[aster_b13_no_data_locations] = self.no_data_value
                ls_emis_data[aster_b14_gap_locations] = 0
                ls_emis_data[aster_b14_no_data_locations] = self.no_data_value

                # Memory cleanup
                del aster_b13_data
                del aster_b14_data
                del aster_b13_gap_locations
                del aster_b14_gap_locations
                del aster_b13_no_data_locations
                del aster_b14_no_data_locations

                # ------------------------------------------------------------
                # Save the no data locations.
                aster_ndvi_no_data_locations = (
                    np.where(aster_ndvi_data == self.no_data_value))

                # Scale the data
                aster_ndvi_data = aster_ndvi_data * 0.01

                # Re-apply the no data locations.
                aster_ndvi_data[aster_ndvi_no_data_locations] = (
                    self.no_data_value)

                # ------------------------------------------------------------
                # Create the estimated Landsat EMIS raster output tile
                try:
                    logger.info('Creating an estimated Landsat EMIS tile')
                    util.Geo.generate_raster_file(gdal.GetDriverByName('GTiff'),
                                                  ls_emis_tile_name,
                                                  ls_emis_data, x_dim, y_dim,
                                                  geo_transform,
                                                  ds_srs.ExportToWkt(),
                                                  self.no_data_value,
                                                  gdal.GDT_Float32)
                except Exception:
                    logger.exception('Creating Landsat EMIS tile')
                    raise

                # ------------------------------------------------------------
                # Create the ASTER NDVI raster output tile
                try:
                    logger.info('Creating an ASTER NDVI tile')
                    util.Geo.generate_raster_file(gdal.GetDriverByName('GTiff'),
                                                  aster_ndvi_tile_name,
                                                  aster_ndvi_data,
                                                  x_dim, y_dim,
                                                  geo_transform,
                                                  ds_srs.ExportToWkt(),
                                                  self.no_data_value,
                                                  gdal.GDT_Float32)
                except Exception:
                    logger.exception('Creating ASTER NDVI tile')
                    raise

                # Memory cleanup
                del ls_emis_data
                del aster_ndvi_data

        # Check to see that we downloaded at least one ASTER tile for
        # processing.
        if len(ls_emis_mean_filenames) == 0:
            raise Exception('No ASTER tiles were downloaded')

        # Save the source proj4 string to use during warping
        src_proj4 = ds_srs.ExportToProj4()

        # Define the temporary names
        ls_emis_mosaic_name = 'landsat_emis_mosaic.tif'
        ls_emis_warped_name = 'landsat_emis_warped.tif'
        aster_ndvi_mosaic_name = 'aster_ndvi_mosaic.tif'
        aster_ndvi_warped_name = 'aster_ndvi_warped.tif'

        # Mosaic the estimated Landsat EMIS tiles into the temp EMIS
        try:
            logger.info('Building mosaic for estimated Landsat EMIS')
            util.Geo.mosaic_tiles_into_one_raster(ls_emis_mean_filenames,
                                                  ls_emis_mosaic_name,
                                                  self.no_data_value)
        except Exception:
            logger.exception('Mosaicing EMIS tiles')
            raise

        # Mosaic the ASTER NDVI tiles into the temp NDVI
        try:
            logger.info('Building mosaic for ASTER NDVI')
            util.Geo.mosaic_tiles_into_one_raster(aster_ndvi_mean_filenames,
                                                  aster_ndvi_mosaic_name,
                                                  self.no_data_value)
        except Exception:
            logger.exception('Mosaicing ASTER NDVI tiles')
            raise

        if not self.keep_intermediate_data:
            # Cleanup the estimated Landsat EMIS tiles
            for emis_filename in ls_emis_mean_filenames:
                if os.path.exists(emis_filename):
                    os.unlink(emis_filename)

            # Cleanup the ASTER NDVI tiles
            for ndvi_filename in aster_ndvi_mean_filenames:
                if os.path.exists(ndvi_filename):
                    os.unlink(ndvi_filename)

        # Warp estimated Landsat EMIS to match the Landsat data
        try:
            logger.info('Warping estimated Landsat EMIS to'
                             ' match Landsat data')
            warp_raster(src_info, src_proj4, self.no_data_value,
                        ls_emis_mosaic_name, ls_emis_warped_name)
        except Exception:
            logger.exception('Warping EMIS to match Landsat data')
            raise

        # Warp ASTER NDVI to match the Landsat data
        try:
            logger.info('Warping ASTER NDVI to match Landsat data')
            warp_raster(src_info, src_proj4, self.no_data_value,
                        aster_ndvi_mosaic_name, aster_ndvi_warped_name)
        except Exception:
            logger.exception('Warping ASTER NDVI to match Landsat data')
            raise

        if not self.keep_intermediate_data:
            # Cleanup the temp files
            if os.path.exists(ls_emis_mosaic_name):
                os.unlink(ls_emis_mosaic_name)
            if os.path.exists(aster_ndvi_mosaic_name):
                os.unlink(aster_ndvi_mosaic_name)

        # Load the warped estimated Landsat EMIS into memory
        ls_emis_data = extract_raster_data(ls_emis_warped_name, 1)
        ls_emis_gap_locations = np.where(ls_emis_data == 0)
        ls_emis_no_data_locations = np.where(ls_emis_data == self.no_data_value)

        # Load the warped ASTER NDVI into memory
        aster_ndvi_data = extract_raster_data(aster_ndvi_warped_name, 1)
        aster_ndvi_gap_locations = np.where(aster_ndvi_data == 0)
        aster_ndvi_no_data_locations = (
            np.where(aster_ndvi_data == self.no_data_value))

        # Turn all negative values to zero
        # Use a realy small value so that we don't have negative zero (-0.0)
        aster_ndvi_data[aster_ndvi_data < 0.0000001] = 0

        if not self.keep_intermediate_data:
            # Cleanup the intermediate files since we have them in memory
            if os.path.exists(ls_emis_warped_name):
                os.unlink(ls_emis_warped_name)
            if os.path.exists(aster_ndvi_warped_name):
                os.unlink(aster_ndvi_warped_name)

        return (ls_emis_data, ls_emis_gap_locations, ls_emis_no_data_locations,
                aster_ndvi_data, aster_ndvi_gap_locations,
                aster_ndvi_no_data_locations)

    # ------------------------------------------------------------------------
    # TODO - NEED TO PROCESS A COASTAL SCENE BECAUSE MAY NOT HAVE ASTER TILE
    #        DATA FOR SOME OF THE SCENE AND WILL NEED TO DO SOMETHING ELSE
    #      - I have implemented a solution for this, but have not been able to
    #        test it.  I have not been able to find a scene that has this
    #        condition.
    # ------------------------------------------------------------------------
    def generate_product(self, espa_metadata):
        '''
        Description:
            Provides the main processing algorithm for generating the
            estimated Landsat emissivity product.  It produces the final
            emissivity product.
        '''

        logger = logging.getLogger(__name__)

        src_info = retrieve_metadata_information(espa_metadata)

        # Determine output information
        dataset = gdal.Open(src_info.toa.red.name)
        output_srs = osr.SpatialReference()
        output_srs.ImportFromWkt(dataset.GetProjection())
        output_transform = dataset.GetGeoTransform()
        output_raster_x_size = dataset.RasterXSize
        output_raster_y_size = dataset.RasterYSize
        del dataset

        coefficients = sensor_coefficients(espa_metadata.xml_object
                                           .global_metadata.satellite)

        # ====================================================================
        # Build NDVI in memory
        (ls_ndvi_data, ndvi_no_data_locations) = (
            generate_landsat_ndvi(src_info, self.no_data_value))

        if self.keep_intermediate_data:
            logger.info('Writing Landsat NDVI raster')
            util.Geo.generate_raster_file(gdal.GetDriverByName('GTiff'),
                                          'internal_landsat_ndvi.tif',
                                          ls_ndvi_data,
                                          output_raster_x_size,
                                          output_raster_y_size,
                                          output_transform,
                                          output_srs.ExportToWkt(),
                                          self.no_data_value,
                                          gdal.GDT_Float32)

        # ====================================================================
        # Determine NDSI and Snow locations
        (snow_locations, ndsi_no_data_locations) = (
            snow_and_ndsi_locations(src_info, self.no_data_value))

        # Build the estimated Landsat EMIS data from the ASTER GED data and
        # warp it to the Landsat scenes projection and image extents
        # For convenience the ASTER NDVI is also extracted and warped to the
        # Landsat scenes projection and image extents
        logger.info('Build thermal emissivity band and retrieve ASTER NDVI')
        (ls_emis_data, ls_emis_gap_locations, ls_emis_no_data_locations,
         aster_ndvi_data, aster_ndvi_gap_locations,
         aster_ndvi_no_data_locations) = self.build_ls_emis_data(src_info,
                                                                 coefficients)

        logger.info('Normalizing Landsat and ASTER NDVI')
        # Normalize Landsat NDVI by max value
        max_ls_ndvi = ls_ndvi_data.max()
        logger.info('Max LS NDVI {0}'.format(max_ls_ndvi))
        ls_ndvi_data = ls_ndvi_data / float(max_ls_ndvi)

        if self.keep_intermediate_data:
            logger.info('Writing Landsat NDVI NORM MAX raster')
            util.Geo.generate_raster_file(gdal.GetDriverByName('GTiff'),
                                          'internal_landsat_ndvi_norm_max.tif',
                                          ls_ndvi_data,
                                          output_raster_x_size,
                                          output_raster_y_size,
                                          output_transform,
                                          output_srs.ExportToWkt(),
                                          self.no_data_value,
                                          gdal.GDT_Float32)

        # Normalize ASTER NDVI by max value
        max_aster_ndvi = aster_ndvi_data.max()
        logger.info('Max ASTER NDVI {0}'.format(max_aster_ndvi))
        aster_ndvi_data = aster_ndvi_data / float(max_aster_ndvi)

        if self.keep_intermediate_data:
            logger.info('Writing Aster NDVI NORM MAX raster')
            util.Geo.generate_raster_file(gdal.GetDriverByName('GTiff'),
                                          'internal_aster_ndvi_norm_max.tif',
                                          aster_ndvi_data,
                                          output_raster_x_size,
                                          output_raster_y_size,
                                          output_transform,
                                          output_srs.ExportToWkt(),
                                          self.no_data_value,
                                          gdal.GDT_Float32)

        # Soil - From prototype code variable name
        logger.info('Calculating EMIS Final')
        with np.errstate(divide='ignore'):
            ls_emis_final = ((ls_emis_data - 0.975 * aster_ndvi_data) /
                             (1.0 - aster_ndvi_data))

        # Memory cleanup
        del aster_ndvi_data
        del ls_emis_data

        # Adjust estimated Landsat EMIS for vegetation and snow, to generate
        # the final Landsat EMIS data
        logger.info('Adjusting estimated EMIS for vegetation')
        ls_emis_final = (coefficients.vegetation_coeff * ls_ndvi_data +
                         ls_emis_final * (1.0 - ls_ndvi_data))

        # Medium snow
        logger.info('Adjusting estimated EMIS for snow')
        ls_emis_final[snow_locations] = coefficients.snow_emissivity

        # Memory cleanup
        del ls_ndvi_data
        del snow_locations

        # Add the fill and scan gaps and ASTER gaps back into the results,
        # since they may have been lost
        logger.info('Adding fill and data gaps back into the estimated'
                         ' Landsat emissivity results')
        ls_emis_final[ls_emis_no_data_locations] = self.no_data_value
        ls_emis_final[ls_emis_gap_locations] = self.no_data_value
        ls_emis_final[aster_ndvi_no_data_locations] = self.no_data_value
        ls_emis_final[aster_ndvi_gap_locations] = self.no_data_value
        ls_emis_final[ndvi_no_data_locations] = self.no_data_value
        ls_emis_final[ndsi_no_data_locations] = self.no_data_value

        # Memory cleanup
        del ls_emis_no_data_locations
        del ls_emis_gap_locations
        del aster_ndvi_no_data_locations
        del aster_ndvi_gap_locations

        product_id = self.xml_filename.split('.xml')[0]
        ls_emis_img_filename = ''.join([product_id, '_emis', '.img'])
        ls_emis_hdr_filename = ''.join([product_id, '_emis', '.hdr'])
        ls_emis_aux_filename = ''.join([ls_emis_img_filename, '.aux', '.xml'])

        logger.info('Creating {0}'.format(ls_emis_img_filename))
        util.Geo.generate_raster_file(gdal.GetDriverByName('ENVI'),
                                      ls_emis_img_filename,
                                      ls_emis_final,
                                      output_raster_x_size,
                                      output_raster_y_size,
                                      output_transform,
                                      output_srs.ExportToWkt(),
                                      self.no_data_value, gdal.GDT_Float32)

        logger.info('Updating {0}'.format(ls_emis_hdr_filename))
        util.Geo.update_envi_header(ls_emis_hdr_filename, self.no_data_value)

        # Remove the *.aux.xml file generated by GDAL
        if os.path.exists(ls_emis_aux_filename):
            os.unlink(ls_emis_aux_filename)

        logger.info('Adding {0} to {1}'.format(ls_emis_img_filename,
                                                    self.xml_filename))
        # Add the estimated Landsat emissivity to the metadata XML
        espa_metadata = Metadata()
        espa_metadata.parse(self.xml_filename)

        # Create an element maker
        em = objectify.ElementMaker(annotate=False,
                                    namespace=None,
                                    nsmap=None)

        sensor_code = product_id[0:3]
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

        emis_band = em.band()
        emis_band.set('product', 'lst_temp')
        emis_band.set('source', source_product)
        emis_band.set('name', 'landsat_emis')
        emis_band.set('category', 'image')
        emis_band.set('data_type', 'FLOAT32')
        emis_band.set('nlines', base_band.attrib['nlines'])
        emis_band.set('nsamps', base_band.attrib['nsamps'])
        emis_band.set('fill_value', str(self.no_data_value))


        emis_band.short_name = em.element('{0}EMIS'.format(sensor_code))
        emis_band.long_name = em.element('Landsat emissivity estimated'
                                         ' from ASTER GED data')
        emis_band.file_name = em.element(ls_emis_img_filename)

        emis_band.pixel_size = base_band.pixel_size

        emis_band.resample_method = em.element('none')
        emis_band.data_units = em.element('Emissivity Coefficient')

        emis_band.valid_range = em.element()
        # TODO TODO TODO - change these to 0.0 and 1.0 after metadata fix
        emis_band.valid_range.set('min', '0')
        emis_band.valid_range.set('max', '1')

        emis_band.app_version = em.element(util.Version.app_version())

        # Get the production date and time in string format
        # Strip the microseconds and add a Z
        date_now = ('{0}Z'.format(datetime.datetime.now()
                                  .strftime('%Y-%m-%dT%H:%M:%S')))
        emis_band.production_date = em.element(date_now)

        # Append the band to the XML
        espa_metadata.xml_object.bands.append(emis_band)

        # Validate the XML
        espa_metadata.validate()

        # Write it to the XML file
        espa_metadata.write(xml_filename=self.xml_filename)


def retrieve_command_line_arguments():
    # Build the command line argument parser
    description = ('Estimates Landsat Emissivity from ASTER GED data')
    parser = ArgumentParser(description=description)

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

    parser.add_argument('--keep-intermediate-data',
                        action='store_true', dest='keep_intermediate_data',
                        required=False, default=False,
                        help='Keep any intermediate data generated')

    parser.add_argument('--debug',
                        action='store_true', dest='debug',
                        required=False, default=False,
                        help='Turn debug messaging on')

    parser.add_argument('--version',
                        action='store_true', dest='version',
                        required=False, default=False,
                        help='Reports the version of the software')

    args = parser.parse_args()

    # Report the version and exit
    if args.version:
        print(util.Version.version_text())
        sys.exit(0)  # EXIT SUCCESS

    # Verify that the --xml parameter was specified
    if args.xml_filename is None:
        raise Exception('--xml must be specified on the command line')

    return args


def main():
    """Generate Landsat EMIS and ASTER NDVI from ASTER GED tiles for the
       specified Landsat scene.
    """

    args = retrieve_command_line_arguments()

    # Check logging level
    debug_level = logging.INFO
    if args.version:
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
        # XML metadata
        espa_metadata = Metadata()
        espa_metadata.parse(args.xml_filename)

        # Register all the gdal drivers
        gdal.AllRegister()

        # Create the processor object
        processor = EstimateLandsatEmissivity(args.xml_filename,
                                              args.aster_ged_server_name,
                                              args.aster_ged_server_path,
                                              args.keep_intermediate_data)

        # Call the main processing routine
        processor.generate_product(espa_metadata)
    except Exception:
        logger.exception('Processing failed')
        sys.exit(1)  # EXIT FAILURE

    logger.info('*** Generate Estimated Landsat Emissivity - Complete ***')


if __name__ == '__main__':
    main()
