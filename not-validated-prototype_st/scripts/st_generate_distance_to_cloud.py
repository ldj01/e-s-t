#! /usr/bin/env python

'''
    File: st_generate_distance_to_cloud.py

    Purpose: Builds a band for distance to cloud based on pixel_qa input.

    Project: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    License: NASA Open Source Agreement 1.3
'''

import os
import sys
import logging
import datetime
from argparse import ArgumentParser
from collections import namedtuple
import numpy as np
from scipy import ndimage
from osgeo import gdal, osr
from lxml import objectify as objectify

import st_utilities as util
from st_exceptions import MissingBandError
from espa import Metadata

PQA_FILL = 0
PQA_CLOUD = 5
PQA_SINGLE_BIT = 0x01             # 00000001


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


def retrieve_command_line_arguments():
    """Read arguments from the command line

    Returns:
        args <arguments>: The arguments read from the command line
    """

    parser = ArgumentParser(description='Builds distance to cloud band')

    parser.add_argument('--version',
                        action='version',
                        version=util.Version.version_text())

    parser.add_argument('--xml',
                        action='store', dest='xml_filename',
                        required=True, default=None,
                        help='The XML metadata file to use')

    parser.add_argument('--debug',
                        action='store_true', dest='debug',
                        required=False, default=False,
                        help='Output debug messages and/or keep debug data')

    args = parser.parse_args()

    # Verify that the --xml parameter was specified
    if args.xml_filename is None:
        raise Exception('--xml must be specified on the command line')

    return args


SourceInfo = namedtuple('SourceInfo', ('proj4', 'qa_filename'))


def retrieve_metadata_information(espa_metadata):
    """Reads required information from the metadata XML file

    Args:
        espa_metadata <espa.Metadata>: XML metadata

    Returns:
        <SourceInfo>: Populated with source information
    """

    pixel_qa_filename = None

    # Find the pixel_qa band to extract information from
    for band in espa_metadata.xml_object.bands.band:
        if (band.get('product') == 'level2_qa' and
                band.get('name') == 'pixel_qa'):
            pixel_qa_filename = str(band.file_name)

            # Get the output proj4 string
            proj4 = util.Geo.get_proj4_projection_string(pixel_qa_filename)


    # Error if we didn't find the required pixel_qa band in the data
    if pixel_qa_filename is None:
        raise MissingBandError('Failed to find the PIXEL QA band'
                               ' in the input data')

    return SourceInfo(proj4=proj4, qa_filename=pixel_qa_filename)


def calculate_distance(src_info, fill_value):
    """Calculate distance to cloud

    Args:
        src_info <SourceInfo>: Information about the source data
        fill_value <float>: No data (fill) value to use

    Returns:
        <numpy.2darray>: Generated distance to cloud band data
    """

    logger = logging.getLogger(__name__)

    logger.info('Building distance to cloud band')

    # Read the pixel QA input
    qa_data = extract_raster_data(src_info.qa_filename, 1)

    # Make a layer that is only set for cloud locations
    qa_cloud_shifted_mask = np.right_shift(qa_data, PQA_CLOUD)
    qa_cloud_mask = np.bitwise_and(qa_cloud_shifted_mask, PQA_SINGLE_BIT)

    # Make cloud pixels the background, and non-cloud the foreground, since
    # distance_transform_edt finds the distance to the closest background pixel
    qa_cloud_background = np.logical_not(qa_cloud_mask)

    # Calculate the distance to clouds
    distance_data = ndimage.distance_transform_edt(qa_cloud_background)

    # Multiply by pixel size in km
    distance_data = distance_data * 0.03

    # Cleanup no data locations.  Skip the right shift since PQA_FILL is 0
    qa_fill_mask = np.bitwise_and(qa_data, PQA_SINGLE_BIT)
    qa_fill_locations = np.where(qa_fill_mask == 1)
    distance_data[qa_fill_locations] = fill_value

    # Memory cleanup
    del qa_data
    del qa_fill_mask
    del qa_fill_locations
    del qa_cloud_shifted_mask
    del qa_cloud_mask
    del qa_cloud_background

    return distance_data


def write_distance_to_cloud_product(samps, lines, transform, wkt, no_data_value,
                                    filename, file_data):
    """Creates the distance to cloud band file

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


def add_cloud_distance_band_to_xml(espa_metadata, filename, sensor_code,
                                   no_data_value):
    """Adds the cloud distance band to the Metadata XML file

    Args:
        espa_metadata <espa.Metadata>: XML metadata information
        filename <str>: Full path for the output file to create
        sensor_code <str>: Name prefix for the sensor
        no_data_value <float>: Value to use for fill
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

    distance_band = maker.band()
    distance_band.set('product', 'st_intermediate')
    distance_band.set('source', source_product)
    distance_band.set('name', 'st_cloud_distance')
    distance_band.set('category', 'image')
    distance_band.set('data_type', 'FLOAT32')
    distance_band.set('nlines', base_band.attrib['nlines'])
    distance_band.set('nsamps', base_band.attrib['nsamps'])
    distance_band.set('fill_value', str(no_data_value))

    distance_band.short_name \
        = maker.element('{0}ST_CLOUD_DIST'.format(sensor_code))

    distance_band.long_name = maker.element('Surface temperature distance to'
                                            ' cloud band')
    distance_band.file_name = maker.element(filename)

    distance_band.pixel_size = base_band.pixel_size

    distance_band.resample_method = maker.element('none')
    distance_band.data_units = maker.element('distance (km)')

    distance_band.valid_range = maker.element()
    distance_band.valid_range.set('min', '0')
    # The largest distance should be about 8000 pixels * 30 meters, or 240 km
    distance_band.valid_range.set('max', '240')

    distance_band.app_version = maker.element(util.Version.app_version())

    # Get the production date and time in string format
    # Strip the microseconds and add a Z
    date_now = ('{0}Z'.format(datetime.datetime.now()
                              .strftime('%Y-%m-%dT%H:%M:%S')))
    distance_band.production_date = maker.element(date_now)

    # Append the band to the XML
    espa_metadata.xml_object.bands.append(distance_band)

    # Validate the XML
    espa_metadata.validate()

    # Write it to the XML file
    espa_metadata.write()


def generate_distance(xml_filename, no_data_value):
    """Provides the main processing algorithm for generating the distance
       to cloud product.

    Args:
        xml_filename <str>: Filename for the ESPA Metadata XML
        no_data_value <float>: No data (fill) value to use
    """

    logger = logging.getLogger(__name__)

    # XML metadata
    espa_metadata = Metadata(xml_filename)
    espa_metadata.parse()

    src_info = retrieve_metadata_information(espa_metadata)

    # Determine output information
    sensor_code = get_satellite_sensor_code(xml_filename)
    dataset = gdal.Open(src_info.qa_filename)
    output_srs = osr.SpatialReference()
    output_srs.ImportFromWkt(dataset.GetProjection())
    output_transform = dataset.GetGeoTransform()
    samps = dataset.RasterXSize
    lines = dataset.RasterYSize
    del dataset

    # Build distance to cloud information in memory
    distance_to_cloud = calculate_distance(src_info, no_data_value)

    # Build distance to cloud filename
    distance_img_filename = ''.join([xml_filename.split('.xml')[0],
                                     '_st_cloud_distance', '.img'])

    # Write distance to cloud product
    write_distance_to_cloud_product(samps=samps,
                                    lines=lines,
                                    transform=output_transform,
                                    wkt=output_srs.ExportToWkt(),
                                    no_data_value=no_data_value,
                                    filename=distance_img_filename,
                                    file_data=distance_to_cloud)

    add_cloud_distance_band_to_xml(espa_metadata=espa_metadata,
                                   filename=distance_img_filename,
                                   sensor_code=sensor_code,
                                   no_data_value=no_data_value)


def get_satellite_sensor_code(xml_filename):
    """Derives the satellite-sensor code from the XML filename

    Args:
        xml_filename <str>: Filename for the ESPA Metadata XML

    Returns:
        <str>: Satellite sensor code
    """

    collection_prefixes = ['LT04', 'LT05', 'LE07', 'LT08', 'LC08', 'LO08']

    base_name = os.path.basename(xml_filename)

    satellite_sensor_code = base_name[0:4]
    if satellite_sensor_code in collection_prefixes:
        return satellite_sensor_code

    raise Exception('Satellite-Sensor code ({0}) not understood'
                    .format(satellite_sensor_code))


# Specify the no data value we will be using.
NO_DATA_VALUE = -9999


def main():
    """Main processing for building the distance to cloud band
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

    logger.info('*** Begin ST Generate Distance to Cloud ***')

    try:
        # Register all the gdal drivers
        gdal.AllRegister()

        # Call the main processing routine
        generate_distance(xml_filename=args.xml_filename,
                          no_data_value=NO_DATA_VALUE)

    except Exception:
        logger.exception('Processing failed')
        sys.exit(1)  # EXIT FAILURE

    logger.info('*** ST Generate Distance to Cloud - Complete ***')

if __name__ == '__main__':
    main()
