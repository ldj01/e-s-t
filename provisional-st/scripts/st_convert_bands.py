#! /usr/bin/env python

'''
    File: st_convert_bands.py

    Purpose: Converts several intermediate bands from float32 (used in science 
             calculations) to int16 to reduce data storage and download 
             requirements

    Project: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    License: NASA Open Source Agreement 1.3

'''

import os
import sys
import logging
from argparse import ArgumentParser
from collections import namedtuple

from osgeo import gdal, osr

import st_utilities as util
from st_exceptions import MissingBandError
from espa import Metadata

SourceInfo = namedtuple('SourceInfo', ('proj4', 'filename'))

MAX_INT16 = 32767   # Maximum GDAL signed integer value 
MIN_INT16 = -32768  # Minimum GDAL signed integer value 

EMIS_SCALE_FACTOR = 0.0001
EMIS_MULT_FACTOR = 10000.0
EMIS_RANGE_MIN = 0.0
EMIS_RANGE_MAX = 10000.0
EMIS_SOURCE_PRODUCT = 'st_intermediate'
EMIS_BAND_NAME = 'emis'

EMIS_STDEV_SCALE_FACTOR = 0.0001
EMIS_STDEV_MULT_FACTOR = 10000.0
EMIS_STDEV_RANGE_MIN = 0.0
EMIS_STDEV_RANGE_MAX = 10000.0
EMIS_STDEV_SOURCE_PRODUCT = 'st_intermediate'
EMIS_STDEV_BAND_NAME = 'emis_stdev'

CLOUD_DISTANCE_SCALE_FACTOR = 0.01
CLOUD_DISTANCE_MULT_FACTOR = 100.0
CLOUD_DISTANCE_RANGE_MIN = 0.0
CLOUD_DISTANCE_RANGE_MAX = 24000.0
CLOUD_DISTANCE_SOURCE_PRODUCT = 'st_intermediate'
CLOUD_DISTANCE_BAND_NAME = 'st_cloud_distance'

# The range is based on the maximum value of the thermal band range with
# gain and bias applied.  The highest value was 21,989 for L8.
THERMAL_RADIANCE_SCALE_FACTOR = 0.001
THERMAL_RADIANCE_MULT_FACTOR = 1000.0
THERMAL_RADIANCE_RANGE_MIN = 0.0
THERMAL_RADIANCE_RANGE_MAX = 22000.0
THERMAL_RADIANCE_SOURCE_PRODUCT = 'st_intermediate'
THERMAL_RADIANCE_BAND_NAME = 'st_thermal_radiance'

UPWELLED_RADIANCE_SCALE_FACTOR = 0.001
UPWELLED_RADIANCE_MULT_FACTOR = 1000.0
UPWELLED_RADIANCE_RANGE_MIN = 0.0
UPWELLED_RADIANCE_RANGE_MAX = 28000.0
UPWELLED_RADIANCE_SOURCE_PRODUCT = 'st_intermediate'
UPWELLED_RADIANCE_BAND_NAME = 'st_upwelled_radiance'

DOWNWELLED_RADIANCE_SCALE_FACTOR = 0.001
DOWNWELLED_RADIANCE_MULT_FACTOR = 1000.0
DOWNWELLED_RADIANCE_RANGE_MIN = 0.0
DOWNWELLED_RADIANCE_RANGE_MAX = 28000.0
DOWNWELLED_RADIANCE_SOURCE_PRODUCT = 'st_intermediate'
DOWNWELLED_RADIANCE_BAND_NAME = 'st_downwelled_radiance'

ATMOSPHERIC_TRANSMITTANCE_SCALE_FACTOR = 0.0001
ATMOSPHERIC_TRANSMITTANCE_MULT_FACTOR = 10000.0
ATMOSPHERIC_TRANSMITTANCE_RANGE_MIN = 0.0
ATMOSPHERIC_TRANSMITTANCE_RANGE_MAX = 10000.0
ATMOSPHERIC_TRANSMITTANCE_SOURCE_PRODUCT = 'st_intermediate'
ATMOSPHERIC_TRANSMITTANCE_BAND_NAME = 'st_atmospheric_transmittance'


def retrieve_command_line_arguments():
    """Read arguments from the command line

    Returns:
        args <arguments>: The arguments read from the command line
    """

    parser = ArgumentParser(description='Convert intermediate bands to int16')

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


def retrieve_metadata_information(espa_metadata, band_name, source_product):
    """Reads required information from the metadata XML file

    Args:
        espa_metadata <espa.Metadata>: XML metadata
        band_name <str>: Name of band to extract metadata information from
        source_product <str>: Source product string for band in the metadata 

    Returns:
        <SourceInfo>: Populated with source information
    """

    intermediate_filename = None

    # Find the intermediate band to extract information from
    for band in espa_metadata.xml_object.bands.band:
        if (band.get('product') == source_product and
                band.get('name') == band_name):
            intermediate_filename = str(band.file_name)

            # Get the output proj4 string
            proj4 = util.Geo.get_proj4_projection_string(intermediate_filename)

    # Error if we didn't find the required intermediate band in the data
    if intermediate_filename is None:
        raise MissingBandError('Failed to find the intermediate band'
                               ' in the input data')

    return SourceInfo(proj4=proj4, filename=intermediate_filename)


def update_band_xml(espa_metadata, source_product, band_name, filename, 
                    scale_factor, range_min, range_max):
    """Updates the intermediate band in the Metadata XML file

    Args:
        espa_metadata <espa.Metadata>: XML metadata information
        source_product <str>: Source product in the metadata 
        band_name <str>: Band name in the metadata 
        filename <str>: Full path for the output file to create
        scale_factor <str>: Value to use for scale factor 
        range_min <str>: Minimum of the data range after scaling 
        range_max <str>: Maximum of the data range after scaling 
    """

    logger = logging.getLogger(__name__)

    logger.info('Updating {0} in Metadata XML'.format(filename))

    # Find the band to update. 
    base_band = None
    for band in espa_metadata.xml_object.bands.band:
        if (band.get('product') == source_product and
                band.get('name') == band_name):
            base_band = band

            band.set('data_type', 'INT16')
            band.set('scale_factor', scale_factor)
            band.valid_range.set('min', range_min)
            band.valid_range.set('max', range_max)

    if base_band is None:
        raise MissingBandError('Failed to find the band in the input data')

    # Validate the XML
    espa_metadata.validate()

    # Write it to the XML file
    espa_metadata.write()


def write_product(samps, lines, transform, wkt, no_data_value, filename,
                  file_data, mult_factor):
    """Creates the converted intermediate band file

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

    # Scale the data
    file_data[file_data != no_data_value] *= mult_factor

    # If the result is outside the new range, put back in the range
    file_data[file_data > MAX_INT16] = MAX_INT16 
    file_data[file_data < MIN_INT16] = MIN_INT16 

    logger = logging.getLogger(__name__)

    logger.info('Updating {0}'.format(filename))
    util.Geo.generate_raster_file(gdal.GetDriverByName('ENVI'),
                                  filename,
                                  file_data,
                                  samps,
                                  lines,
                                  transform,
                                  wkt,
                                  no_data_value,
                                  gdal.GDT_Int16)

    hdr_filename = filename.replace('.img', '.hdr')
    logger.info('Updating {0}'.format(hdr_filename))
    util.Geo.update_envi_header(hdr_filename, no_data_value)

    # Remove the *.aux.xml file generated by GDAL
    aux_filename = filename.replace('.img', '.img.aux.xml')
    if os.path.exists(aux_filename):
        os.unlink(aux_filename)


def convert_band(espa_metadata, xml_filename, no_data_value, scale_factor, 
                 mult_factor, range_min, range_max, source_product, band_name):
    """Convert a single intermediate band

    Args:
        espa_metadata <espa.Metadata>: XML metadata
        xml_filename <str>: Filename for the ESPA Metadata XML
        no_data_value <float>: No data (fill) value to use
        scale_factor <str>: Scale factor used in the conversion 
        mult_factor <float>: Multiplication factor used in the conversion 
        range_min <str>: Minimum of the data range after scaling 
        range_max <str>: Maximum of the data range after scaling 
        source_product <str>: Source product string for band in the metadata 
        band_name <str>: Band name string in the metadata 
    """

    # Determine output information.
    src_info = retrieve_metadata_information(espa_metadata, band_name,
                                             source_product)
    dataset = gdal.Open(src_info.filename)
    output_srs = osr.SpatialReference()
    output_srs.ImportFromWkt(dataset.GetProjection())
    output_transform = dataset.GetGeoTransform()
    samps = dataset.RasterXSize
    lines = dataset.RasterYSize
    del dataset

    # Read band
    data_array = util.Dataset.extract_raster_data(src_info.filename, 1)

    # Build converted intermediate band filename
    product_id = espa_metadata.xml_object.global_metadata.product_id.text
    img_filename = ''.join([product_id, '_' + band_name, '_converted.img'])

    # Write updated intermediate product
    write_product(samps=samps,
                  lines=lines,
                  transform=output_transform,
                  wkt=output_srs.ExportToWkt(),
                  no_data_value=no_data_value,
                  filename=img_filename,
                  file_data=data_array,
                  mult_factor=mult_factor)

    # Update the band's metadata to reflect conversion changes
    update_band_xml(espa_metadata=espa_metadata,
                    source_product=source_product,
                    band_name=band_name,
                    filename=img_filename,
                    scale_factor=scale_factor,
                    range_min=range_min,
                    range_max=range_max)


def convert_bands(xml_filename, no_data_value):
    """Convert multiple intermediate bands

    Args:
        xml_filename <str>: Filename for the ESPA Metadata XML
        no_data_value <float>: No data (fill) value to use
    """

    # XML metadata
    espa_metadata = Metadata(xml_filename)
    espa_metadata.parse()

    # Convert emissivity band.
    convert_band(espa_metadata=espa_metadata,
                 xml_filename=xml_filename,
                 no_data_value=no_data_value,
                 scale_factor=str(EMIS_SCALE_FACTOR),
                 mult_factor=EMIS_MULT_FACTOR,
                 range_min=str(EMIS_RANGE_MIN),
                 range_max=str(EMIS_RANGE_MAX),
                 source_product=EMIS_SOURCE_PRODUCT,
                 band_name=EMIS_BAND_NAME)

    # Convert emissivity standard deviation band.
    convert_band(espa_metadata=espa_metadata,
                 xml_filename=xml_filename,
                 no_data_value=no_data_value,
                 scale_factor=str(EMIS_STDEV_SCALE_FACTOR),
                 mult_factor=EMIS_STDEV_MULT_FACTOR,
                 range_min=str(EMIS_STDEV_RANGE_MIN),
                 range_max=str(EMIS_STDEV_RANGE_MAX),
                 source_product=EMIS_STDEV_SOURCE_PRODUCT,
                 band_name=EMIS_STDEV_BAND_NAME)

    # Convert cloud distance band.
    convert_band(espa_metadata=espa_metadata,
                 xml_filename=xml_filename,
                 no_data_value=no_data_value,
                 scale_factor=str(CLOUD_DISTANCE_SCALE_FACTOR),
                 mult_factor=CLOUD_DISTANCE_MULT_FACTOR,
                 range_min=str(CLOUD_DISTANCE_RANGE_MIN),
                 range_max=str(CLOUD_DISTANCE_RANGE_MAX),
                 source_product=CLOUD_DISTANCE_SOURCE_PRODUCT,
                 band_name=CLOUD_DISTANCE_BAND_NAME)

    # Convert thermal radiance band.
    convert_band(espa_metadata=espa_metadata,
                 xml_filename=xml_filename,
                 no_data_value=no_data_value,
                 scale_factor=str(THERMAL_RADIANCE_SCALE_FACTOR),
                 mult_factor=THERMAL_RADIANCE_MULT_FACTOR,
                 range_min=str(THERMAL_RADIANCE_RANGE_MIN),
                 range_max=str(THERMAL_RADIANCE_RANGE_MAX),
                 source_product=THERMAL_RADIANCE_SOURCE_PRODUCT,
                 band_name=THERMAL_RADIANCE_BAND_NAME)

    # Convert upwelled radiance band.
    convert_band(espa_metadata=espa_metadata,
                 xml_filename=xml_filename,
                 no_data_value=no_data_value,
                 scale_factor=str(UPWELLED_RADIANCE_SCALE_FACTOR),
                 mult_factor=UPWELLED_RADIANCE_MULT_FACTOR,
                 range_min=str(UPWELLED_RADIANCE_RANGE_MIN),
                 range_max=str(UPWELLED_RADIANCE_RANGE_MAX),
                 source_product=UPWELLED_RADIANCE_SOURCE_PRODUCT,
                 band_name=UPWELLED_RADIANCE_BAND_NAME)

    # Convert downwelled radiance band.
    convert_band(espa_metadata=espa_metadata,
                 xml_filename=xml_filename,
                 no_data_value=no_data_value,
                 scale_factor=str(DOWNWELLED_RADIANCE_SCALE_FACTOR),
                 mult_factor=DOWNWELLED_RADIANCE_MULT_FACTOR,
                 range_min=str(DOWNWELLED_RADIANCE_RANGE_MIN),
                 range_max=str(DOWNWELLED_RADIANCE_RANGE_MAX),
                 source_product=DOWNWELLED_RADIANCE_SOURCE_PRODUCT,
                 band_name=DOWNWELLED_RADIANCE_BAND_NAME)

    # Convert atmospheric transmittance band.
    convert_band(espa_metadata=espa_metadata,
                 xml_filename=xml_filename,
                 no_data_value=no_data_value,
                 scale_factor=str(ATMOSPHERIC_TRANSMITTANCE_SCALE_FACTOR),
                 mult_factor=ATMOSPHERIC_TRANSMITTANCE_MULT_FACTOR,
                 range_min=str(ATMOSPHERIC_TRANSMITTANCE_RANGE_MIN),
                 range_max=str(ATMOSPHERIC_TRANSMITTANCE_RANGE_MAX),
                 source_product=ATMOSPHERIC_TRANSMITTANCE_SOURCE_PRODUCT,
                 band_name=ATMOSPHERIC_TRANSMITTANCE_BAND_NAME)


def main():
    """Main processing for converting surface temperature intermediate bands
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

    logger.info('*** Begin ST Convert Bands ***')

    try:
        # Register all the gdal drivers
        gdal.AllRegister()

        # Call the main processing routine
        convert_bands(xml_filename=args.xml_filename,
                      no_data_value=util.INTERMEDIATE_NO_DATA_VALUE)

    except Exception:
        logger.exception('Processing failed')
        sys.exit(1)  # EXIT FAILURE

    logger.info('*** ST Convert Bands - Complete ***')

if __name__ == '__main__':
    main()
