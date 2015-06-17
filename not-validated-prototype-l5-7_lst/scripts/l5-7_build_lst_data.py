#! /usr/bin/env python

'''
    FILE: l5-7_build_lst_data.py

    PURPOSE: Builds the LST product from the intermediate data that was
             generated.

    PROJECT: Land Satellites Data Systems (LSDS) Science Research and
             Development (LSRD) at the USGS EROS

    LICENSE: NASA Open Source Agreement 1.3

    Algorithm Authors:

        Monica Cook
        Chester F. Carlson Center for Imaging Science
        Rochester Institute of Technology
        mxc7441@rit.edu

        Glynn C. Hulley
        Research Scientist
        NASA Jet Propulsion Laboratory
        email: glynn.hulley@nasa.gov

        Simon J. Hook
        Research Scientist
        NASA Jet Propulsion Laboratory
        email: Simon.J.Hook@jpl.nasa.gov

    Date              Programmer               Reason
    ----------------  ------------------------ -------------------------------
    June/2015         Ron Dilley               Initial implementation
'''

import os
import sys
import glob
import logging
import math
import requests
import commands
import datetime
import numpy as np
from cStringIO import StringIO
from argparse import ArgumentParser
from osgeo import gdal, osr
from time import sleep


# Import the metadata api found in the espa-product-formatter project
import metadata_api


# Import local modules
import lst_utilities as util


'''
 Specify the no data value we will be using, it also matches the no_data_value
 for the ASTER data we extract and use
'''
NO_DATA_VALUE = -9999


# ============================================================================
def update_envi_header(hdr_file_path, no_data_value):
    '''
    Description:
        Updates the specified ENVI header.  Especially the no data value,
        since it is not supported by the GDAL ENVI driver.
    '''

    sb = StringIO()
    with open(hdr_file_path, 'r') as tmp_fd:
        while True:
            line = tmp_fd.readline()
            if not line:
                break
            if (line.startswith('data ignore value') or
                    line.startswith('description')):
                pass
            else:
                sb.write(line)

            if line.startswith('description'):
                # This may be on multiple lines so read lines until
                # we find the closing brace
                if not line.strip().endswith('}'):
                    while 1:
                        next_line = tmp_fd.readline()
                        if (not next_line or
                                next_line.strip().endswith('}')):
                            break
                sb.write('description = {ESPA-generated file}\n')
            elif (line.startswith('data type') and
                  (no_data_value is not None)):
                sb.write('data ignore value = %s\n' % no_data_value)

    # Do the actual replace here
    with open(hdr_file_path, 'w') as tmp_fd:
        tmp_fd.write(sb.getvalue())


# ============================================================================
def generate_raster_file(driver, filename, x_dim, y_dim, data,
                         geo_transform, proj_wkt, no_data_value, data_type):
    '''
    Description:
        Creates a raster file on disk for the specified data, using the
        specified driver.

    Note: It is assumed that the driver supports setting of the no data value.
          It is the callers responsibility to fix it if it does not.

    Note: It is assumed that the caller specified the correct file extension
          in the filename parameter for the specfied driver.
    '''

    try:
        raster = driver.Create(filename, x_dim, y_dim, 1, data_type)

        raster.SetGeoTransform(geo_transform)
        raster.SetProjection(proj_wkt)
        raster.GetRasterBand(1).WriteArray(data)
        raster.GetRasterBand(1).SetNoDataValue(no_data_value)
        raster.FlushCache()

        # Cleanup memory
        del (raster)

    except Exception:
        raise


# ============================================================================
# OK OK OK
def retrieve_metadata_information(xml_filename):
    '''
    Description:
        Loads and reads required information from the metadata XML file.
    '''

    # Read the XML metadata
    espa_xml = metadata_api.parse(xml_filename, silence=True)
    # Grab the global metadata object
    gm = espa_xml.get_global_metadata()
    # Grab the bands metadata object
    bands = espa_xml.get_bands()

    thermal_name = ''
    transmittance_name = ''
    upwelled_name = ''
    downwelled_name = ''
    emissivity_name = ''

    # Find the TOA bands to extract information from
    for band in bands.band:
        if (band.product == 'lst_temp' and
                band.name == 'lst_thermal_radiance'):
            thermal_name = band.get_file_name()

        if (band.product == 'lst_temp' and
                band.name == 'lst_atmospheric_transmittance'):
            transmittance_name = band.get_file_name()

        if (band.product == 'lst_temp' and
                band.name == 'lst_upwelled_radiance'):
            upwelled_name = band.get_file_name()

        if (band.product == 'lst_temp' and
                band.name == 'lst_downwelled_radiance'):
            downwelled_name = band.get_file_name()

        if (band.product == 'lst_temp' and
                band.name == 'landsat_emis'):
            emissivity_name = band.get_file_name()

    # Error if we didn't find the required TOA bands in the data
    if len(thermal_name) <= 0:
        raise Exception("Failed to find the lst_thermal_radiance band in the"
                        " input data")
    if len(transmittance_name) <= 0:
        raise Exception("Failed to find the lst_atmospheric_transmittance in"
                        " the input data")
    if len(upwelled_name) <= 0:
        raise Exception("Failed to find the lst_upwelled_radiance in the"
                        " input data")
    if len(downwelled_name) <= 0:
        raise Exception("Failed to find the lst_downwelled_radiance in the"
                        " input data")
    if len(emissivity_name) <= 0:
        raise Exception("Failed to find the landsat_emis in the input data")

    # Save for later
    satellite = gm.satellite

    del (bands)
    del (gm)
    del (espa_xml)

    return (satellite, thermal_name, transmittance_name, upwelled_name,
            downwelled_name, emissivity_name)


# ============================================================================
def process(args, lst_data_dir):
    '''
    Description:
        Provides the main processing algorithm for building the Land Surface
        Temperature product.  It produces the final LST product.
    '''

    logger = logging.getLogger(__name__)

    try:
        (satellite, thermal_name, trans_name, upwelled_name, downwelled_name,
         emissivity_name) = retrieve_metadata_information(args.xml_filename)
    except:
        logger.exception("Failed reading input XML metadata file")
        raise

    # Register all the gdal drivers and choose the ENVI for our output
    gdal.AllRegister()
    envi_driver = gdal.GetDriverByName('ENVI')

    # Read the bands into memory
    logger.info("Loading intermediate data")
    ds = gdal.Open(thermal_name)
    x_dim = ds.RasterXSize  # They are all the same size
    y_dim = ds.RasterYSize
    thermal_data = ds.GetRasterBand(1).ReadAsArray(0, 0, x_dim, y_dim)
    thermal_masked = np.ma.masked_equal(thermal_data, NO_DATA_VALUE)

    ds = gdal.Open(trans_name)
    trans_data = ds.GetRasterBand(1).ReadAsArray(0, 0, x_dim, y_dim)
    trans_masked = np.ma.masked_equal(trans_data, NO_DATA_VALUE)

    ds = gdal.Open(upwelled_name)
    upwelled_data = ds.GetRasterBand(1).ReadAsArray(0, 0, x_dim, y_dim)
    upwelled_masked = np.ma.masked_equal(upwelled_data, NO_DATA_VALUE)

    ds = gdal.Open(downwelled_name)
    downwelled_data = ds.GetRasterBand(1).ReadAsArray(0, 0, x_dim, y_dim)
    downwelled_masked = np.ma.masked_equal(downwelled_data, NO_DATA_VALUE)

    ds = gdal.Open(emissivity_name)
    emissivity_data = ds.GetRasterBand(1).ReadAsArray(0, 0, x_dim, y_dim)
    emissivity_masked = np.ma.masked_equal(emissivity_data, NO_DATA_VALUE)

    # Save the locations of the fill and scan gaps
    thermal_no_data_locations = np.where(thermal_data == NO_DATA_VALUE)
    trans_no_data_locations = np.where(trans_data == NO_DATA_VALUE)
    upwelled_no_data_locations = np.where(upwelled_data == NO_DATA_VALUE)
    downwelled_no_data_locations = np.where(downwelled_data == NO_DATA_VALUE)
    emissivity_no_data_locations = np.where(emissivity_data == NO_DATA_VALUE)

    # Save for the output product
    ds_srs = osr.SpatialReference()
    ds_srs.ImportFromWkt(ds.GetProjection())
    ds_transform = ds.GetGeoTransform()

    # Memory cleanup
    del (thermal_data)
    del (trans_data)
    del (upwelled_data)
    del (downwelled_data)
    del (emissivity_data)
    del (ds)

    # Build the LST data from the intermediate inputs
    logger.info("Calculating surface radiance")
    surface_radiance = (thermal_masked - upwelled_masked) / trans_masked
    radiance = (surface_radiance -
                (1.0 - emissivity_masked) * downwelled_masked)

    # Account for emissivity to get Plank emitted radiance
    logger.info("Calculating Plank emitted radiance")
    radiance_emitted = radiance / emissivity_masked

    # Memory cleanup
    del (thermal_masked)
    del (trans_masked)
    del (upwelled_masked)
    del (downwelled_masked)
    del (emissivity_masked)
    del (surface_radiance)
    del (radiance)

    # Use Brightness Temperature LUT to get skin temperature
    # Read the correct one for what we are processing
    if satellite == 'LANDSAT_7':
        logger.info("Using Landsat 7 Brightness Temperature LUT")
        bt_name = "L7_Brightness_Temperature_LUT.txt"

    elif satellite == 'LANDSAT_5':
        logger.info("Using Landsat 5 Brightness Temperature LUT")
        bt_name = "L5_Brightness_Temperature_LUT.txt"

    bt_data = np.loadtxt(os.path.join(lst_data_dir, bt_name),
                         dtype=float, delimiter=' ')
    bt_radiance_LUT = bt_data[:, 1]
    bt_temp_LUT = bt_data[:, 0]

    logger.info("Generating LST results")
    lst_data = np.interp(radiance_emitted, bt_radiance_LUT, bt_temp_LUT)

    # Memory cleanup
    del (radiance_emitted)

    # TODO TODO TODO - SHOULD PROBABLY SCALE/CONVERT THE DATA TO INT16

    # Add the fill and scan gaps back into the results, since they may have
    # been lost
    logger.info("Adding fill and data gaps back into the Land Surface"
                " Temperature results")
    lst_data[thermal_no_data_locations] = NO_DATA_VALUE
    lst_data[trans_no_data_locations] = NO_DATA_VALUE
    lst_data[upwelled_no_data_locations] = NO_DATA_VALUE
    lst_data[downwelled_no_data_locations] = NO_DATA_VALUE
    lst_data[emissivity_no_data_locations] = NO_DATA_VALUE

    # Memory cleanup
    del (thermal_no_data_locations)
    del (trans_no_data_locations)
    del (upwelled_no_data_locations)
    del (downwelled_no_data_locations)
    del (emissivity_no_data_locations)

    product_id = args.xml_filename.split('.xml')[0]
    lst_img_filename = ''.join([product_id, '_lst', '.img'])
    lst_hdr_filename = ''.join([product_id, '_lst', '.hdr'])
    lst_aux_filename = ''.join([lst_img_filename, '.aux', '.xml'])

    logger.info("Creating {0}".format(lst_img_filename))
    generate_raster_file(envi_driver, lst_img_filename, x_dim, y_dim,
                         lst_data, ds_transform,
                         ds_srs.ExportToWkt(), NO_DATA_VALUE,
                         gdal.GDT_Float32)
    # TODO TODO TODO - FIX ABOVE IF CONVERTED TO INT16

    logger.info("Updating {0}".format(lst_hdr_filename))
    update_envi_header(lst_hdr_filename, NO_DATA_VALUE)

    # Remove the *.aux.xml file generated by GDAL
    if os.path.exists(lst_aux_filename):
        os.unlink(lst_aux_filename)

    # Add the estimated Land Surface Temperature product to the metadata XML
    espa_xml = metadata_api.parse(args.xml_filename, silence=True)
    bands = espa_xml.get_bands()
    sensor_code = product_id[0:3]

    # Find the TOA Band 1 to use for the specific band details
    base_band = None
    for band in bands.band:
        if band.product == 'toa_refl' and band.name == 'toa_band1':
            base_band = band

    if base_band is None:
        raise Exception("Failed to find the TOA BLUE band in the input data")

    lst_band = metadata_api.band(product="lst",
                                 source='toa_refl',
                                 name="land_surface_temperature",
                                 category="image",
                                 data_type="FLOAT32",
    # TODO TODO TODO - FIX ABOVE IF CONVERTED TO INT16
                                 nlines=base_band.get_nlines(),
                                 nsamps=base_band.get_nsamps(),
                                 fill_value=str(NO_DATA_VALUE))

    lst_band.set_short_name('{0}LST'.format(sensor_code))
    lst_band.set_long_name("Land Surface Temperature")
    lst_band.set_file_name(lst_img_filename)
    lst_band.set_data_units("temperature (kelvin)")

    pixel_size = metadata_api.pixel_size(base_band.pixel_size.x,
                                         base_band.pixel_size.x,
                                         base_band.pixel_size.units)
    lst_band.set_pixel_size(pixel_size)

    valid_range = metadata_api.valid_range(min=150.0, max=373.0)
    # TODO TODO TODO - FIX ABOVE IF CONVERTED TO INT16
    lst_band.set_valid_range(valid_range)

    # Set the date, but first clean the microseconds off of it
    production_date = \
        datetime.datetime.strptime(datetime.datetime.now().
                                   strftime('%Y-%m-%dT%H:%M:%S'),
                                   '%Y-%m-%dT%H:%M:%S')

    lst_band.set_production_date(production_date)

    lst_band.set_app_version(util.Version.app_version())

    bands.add_band(lst_band)

    # Write the XML metadata file out
    with open(args.xml_filename, 'w') as fd:
        metadata_api.export(fd, espa_xml)

    # Memory cleanup
    del (lst_data)


# ============================================================================
if __name__ == '__main__':
    '''
    Description:
        Generate Landsat EMIS and ASTER NDVI from ASTER GED tiles for the
        specified Landsat scene.
    '''

    # Build the command line argument parser
    description = ("Reads intermediate data generated previously and combines"
                   " them into the Land Surface Temperature product")
    parser = ArgumentParser(description=description)

    # ---- Add parameters ----
    # Required parameters
    parser.add_argument('--xml',
                        action='store', dest='xml_filename',
                        required=False, default=None,
                        help="The XML metadata file to use")

    # Optional parameters
    parser.add_argument('--version',
                        action='store_true', dest='version',
                        required=False, default=False,
                        help="Reports the version of the software")

    # Parse the command line arguments
    args = parser.parse_args()

    # Command line arguments are required so print the help if none were
    # provided
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)  # EXIT FAILURE

    # Report the version and exit
    if args.version:
        print util.Version.version_text()
        sys.exit(0)  # EXIT SUCCESS

    # Verify that the --xml parameter was specified
    if args.xml_filename is None:
        raise Exception("--xml must be specified on the command line")
        sys.exit(1)  # EXIT FAILURE

    # Configure logging
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:%(funcName)s'
                                ' -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)

    logger = logging.getLogger(__name__)

    # Verify environment variables exist
    lst_data_dir = os.environ.get('LST_DATA_DIR')
    if lst_data_dir is None:
        logger.info("Missing environment variable LST_DATA_DIR")
        sys.exit(1)  # EXIT FAILURE

    try:
        # Call the main processing routine
        process(args, lst_data_dir)
    except Exception, e:
        if hasattr(e, 'output'):
            logger.error("Output [%s]" % e.output)
        logger.exception("Processing failed")
        sys.exit(1)  # EXIT FAILURE

    sys.exit(0)  # EXIT SUCCESS
