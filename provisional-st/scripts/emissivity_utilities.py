#! /usr/bin/env python

'''
    FILE: emissivity_utilities.py

    PURPOSE: Provide a library of routines to be used by ST emissivity
             applications.

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
import logging
import math
import datetime
from argparse import ArgumentParser
from collections import namedtuple

import requests
from lxml import objectify as objectify
from osgeo import gdal

from st_exceptions import MissingBandError


# Import local modules
import st_utilities as util


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
                      ('name', 'scale_factor', 'add_offset', 'pixel_size'))
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
                    add_offset=float(band.get('add_offset')),
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

    cmd = ['gdalwarp', '-wm', '2048', '-wo', 'NUM_THREADS=2',
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


def shift_longitude(tile_name, shifted_tile_name, offset):
    """Shift the longitude of the tile data and put the results in the 
       requested output file

    Args:
        tile_name <str>: Filename of tile to shift 
        shifted_tile_name <str>: Filename of output tile to place results 
        offset <int>: Number of degrees to add to the longitude 
    """

    logger = logging.getLogger(__name__)

    # Set up the base command
    cmd = ['gdal_translate', '-a_ullr']

    # Get the current locations
    tile_src = gdal.Open(tile_name)
    ulx, xres, xskew, uly, yskew, yres = tile_src.GetGeoTransform()

    # Compute the adjusted longitude locations
    lrx = ulx + (tile_src.RasterXSize * xres)
    lry = uly + (tile_src.RasterYSize * yres)
    new_ulx = ulx + offset
    new_lrx = lrx + offset

    # Close the dataset
    tile_src = None

    # Add updated coordinates to the command
    cmd.extend([str(new_ulx), str(uly), str(new_lrx), str(lry)])

    # Add source and destination files to the command
    cmd.extend([tile_name, shifted_tile_name])

    # Convert to a string for the execution
    cmd = ' '.join(cmd)

    output = ''
    try:
        logger.info('Executing [{0}]'.format(cmd))
        output = util.System.execute_cmd(cmd)
    finally:
        if len(output) > 0:
            logger.info(output)


# Shift tile longitudes
def shift_tiles(tiles):
    """Shift the longitude of the tiles that need it to ensure they are in
       the 0..360 range.  This is intended to be used in antimeridian crossing
       cases.  Making the longitude range contiguous enables mosaicking using
       GDAL tools.

    Args:
        tiles <list(<str>)>: List of tiles to shift 
    """

    logger = logging.getLogger(__name__)

    for tile in tiles:
        longitude = int(tile.split(".")[3])

        # Only shift the longitude of the tiles with negative longitude
        if longitude < 0:

            # Name the shifted output file
            shifted_tile = tile + '_shifted'

            # Shift the longitude values
            shift_longitude(tile, shifted_tile, 360)

            # Move destination file back to source file
            output = ''
            try:
                cmd = 'mv {0} {1}'.format(shifted_tile, tile)
                logger.info('Executing [{0}]'.format(cmd))
                output = util.System.execute_cmd(cmd)
            finally:
                if len(output) > 0:
                    logger.info(output)


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
    aux_filename = filename.replace('.img', '.img.aux.xml')
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
        raise MissingBandError('Failed to find the TOA band 1'
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


def get_env_var(variable, default):
    """Looks up the requested environment variable

    Args:
        variable <str>: Environment variable to get
        default <str,int,None>: Default value for the environment variable

    Returns:
        <str>: Value of the environment variable
    """

    result = os.environ.get(variable, default)
    if result == None:
        raise RuntimeError(
            'You must specify {} in the environment'.format(variable))

    return result

