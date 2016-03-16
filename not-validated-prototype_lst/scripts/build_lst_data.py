#! /usr/bin/env python

'''
    FILE: build_lst_data.py

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


SCALE_FACTOR = 0.1
MULT_FACTOR = 10.0


# ============================================================================
class BuildLSTData(object):
    '''
    Description:
        Defines the processor for generating the Land Surface Temperature
        product.
    '''

    def __init__(self, xml_filename):
        super(BuildLSTData, self).__init__()

        # Keep local copies of this
        self.xml_filename = xml_filename

        self.lst_data_dir = ''
        # Grab the data directory from the environment
        if 'LST_DATA_DIR' not in os.environ:
            raise Exception('Environment variable LST_DATA_DIR is'
                            ' not defined')
        else:
            self.lst_data_dir = os.environ.get('LST_DATA_DIR')

        # Specify the no data value we will be using, it also matches the
        # no_data_value for the ASTER data we extract and use
        self.no_data_value = -9999

        # Setup the logger to use
        self.logger = logging.getLogger(__name__)

    # ------------------------------------------------------------------------
    def retrieve_metadata_information(self):
        '''
        Description:
            Loads and reads required information from the metadata XML file.
        '''

        # Read the XML metadata
        espa_xml = metadata_api.parse(self.xml_filename, silence=True)
        # Grab the global metadata object
        gm = espa_xml.get_global_metadata()
        # Grab the bands metadata object
        bands = espa_xml.get_bands()

        self.thermal_name = ''
        self.transmittance_name = ''
        self.upwelled_name = ''
        self.downwelled_name = ''
        self.emissivity_name = ''

        # Find the TOA bands to extract information from
        for band in bands.band:
            if (band.product == 'lst_temp' and
                    band.name == 'lst_thermal_radiance'):
                self.thermal_name = band.get_file_name()

            if (band.product == 'lst_temp' and
                    band.name == 'lst_atmospheric_transmittance'):
                self.transmittance_name = band.get_file_name()

            if (band.product == 'lst_temp' and
                    band.name == 'lst_upwelled_radiance'):
                self.upwelled_name = band.get_file_name()

            if (band.product == 'lst_temp' and
                    band.name == 'lst_downwelled_radiance'):
                self.downwelled_name = band.get_file_name()

            if (band.product == 'lst_temp' and
                    band.name == 'landsat_emis'):
                self.emissivity_name = band.get_file_name()

        # Error if we didn't find the required TOA bands in the data
        if len(self.thermal_name) <= 0:
            raise Exception('Failed to find the lst_thermal_radiance band'
                            ' in the input data')
        if len(self.transmittance_name) <= 0:
            raise Exception('Failed to find the lst_atmospheric_transmittance'
                            ' in the input data')
        if len(self.upwelled_name) <= 0:
            raise Exception('Failed to find the lst_upwelled_radiance'
                            ' in the input data')
        if len(self.downwelled_name) <= 0:
            raise Exception('Failed to find the lst_downwelled_radiance'
                            ' in the input data')
        if len(self.emissivity_name) <= 0:
            raise Exception('Failed to find the landsat_emis'
                            ' in the input data')

        # Save for later
        self.satellite = gm.satellite

        del (bands)
        del (gm)
        del (espa_xml)

    # ------------------------------------------------------------------------
    def generate_data(self):
        '''
        Description:
            Provides the main processing algorithm for building the Land
            Surface Temperature product.  It produces the final LST product.
        '''

        try:
            self.retrieve_metadata_information()
        except Exception:
            self.logger.exception('Failed reading input XML metadata file')
            raise

        # Register all the gdal drivers and choose the ENVI for our output
        gdal.AllRegister()
        envi_driver = gdal.GetDriverByName('ENVI')

        # Read the bands into memory

        # Landsat Radiance at sensor for thermal band
        self.logger.info('Loading intermediate thermal band data')
        ds = gdal.Open(self.thermal_name)
        x_dim = ds.RasterXSize  # They are all the same size
        y_dim = ds.RasterYSize
        thermal_data = ds.GetRasterBand(1).ReadAsArray(0, 0, x_dim, y_dim)

        # Atmospheric transmittance
        self.logger.info('Loading intermediate transmittance band data')
        ds = gdal.Open(self.transmittance_name)
        trans_data = ds.GetRasterBand(1).ReadAsArray(0, 0, x_dim, y_dim)

        # Atmospheric path radiance - upwelled radiance
        self.logger.info('Loading intermediate upwelled band data')
        ds = gdal.Open(self.upwelled_name)
        upwelled_data = ds.GetRasterBand(1).ReadAsArray(0, 0, x_dim, y_dim)

        self.logger.info('Calculating surface radiance')
        # Surface radiance
        surface_radiance = (thermal_data - upwelled_data) / trans_data

        # Fix the no data locations
        no_data_locations = np.where(thermal_data == self.no_data_value)
        surface_radiance[no_data_locations] = self.no_data_value

        no_data_locations = np.where(trans_data == self.no_data_value)
        surface_radiance[no_data_locations] = self.no_data_value

        no_data_locations = np.where(upwelled_data == self.no_data_value)
        surface_radiance[no_data_locations] = self.no_data_value

        # Memory cleanup
        del (thermal_data)
        del (trans_data)
        del (upwelled_data)
        del (no_data_locations)

        # Downwelling sky irradiance
        self.logger.info('Loading intermediate downwelled band data')
        ds = gdal.Open(self.downwelled_name)
        downwelled_data = ds.GetRasterBand(1).ReadAsArray(0, 0, x_dim, y_dim)

        # Landsat emissivity estimated from ASTER GED data
        self.logger.info('Loading intermediate emissivity band data')
        ds = gdal.Open(self.emissivity_name)
        emissivity_data = ds.GetRasterBand(1).ReadAsArray(0, 0, x_dim, y_dim)

        # Save for the output product
        ds_srs = osr.SpatialReference()
        ds_srs.ImportFromWkt(ds.GetProjection())
        ds_transform = ds.GetGeoTransform()

        # Memory cleanup
        del (ds)

        # Estimate Earth-emitted radiance by subtracting off the reflected
        # downwelling component
        radiance = (surface_radiance -
                    (1.0 - emissivity_data) * downwelled_data)

        # Account for surface emissivity to get Plank emitted radiance
        self.logger.info('Calculating Plank emitted radiance')
        radiance_emitted = radiance / emissivity_data

        # Fix the no data locations
        no_data_locations = np.where(surface_radiance == self.no_data_value)
        radiance_emitted[no_data_locations] = self.no_data_value

        no_data_locations = np.where(downwelled_data == self.no_data_value)
        radiance_emitted[no_data_locations] = self.no_data_value

        no_data_locations = np.where(emissivity_data == self.no_data_value)
        radiance_emitted[no_data_locations] = self.no_data_value

        # Memory cleanup
        del (downwelled_data)
        del (emissivity_data)
        del (surface_radiance)
        del (radiance)
        del (no_data_locations)

        # Use Brightness Temperature LUT to get skin temperature
        # Read the correct one for what we are processing
        if self.satellite == 'LANDSAT_7':
            self.logger.info('Using Landsat 7 Brightness Temperature LUT')
            bt_name = 'L7_Brightness_Temperature_LUT.txt'

        elif self.satellite == 'LANDSAT_5':
            self.logger.info('Using Landsat 5 Brightness Temperature LUT')
            bt_name = 'L5_Brightness_Temperature_LUT.txt'

        bt_data = np.loadtxt(os.path.join(self.lst_data_dir, bt_name),
                             dtype=float, delimiter=' ')
        bt_radiance_LUT = bt_data[:, 1]
        bt_temp_LUT = bt_data[:, 0]

        self.logger.info('Generating LST results')
        lst_data = np.interp(radiance_emitted, bt_radiance_LUT, bt_temp_LUT)

        # Scale the result and convert it to an int16
        lst_data = lst_data * MULT_FACTOR
        lst_fata = lst_data.astype(np.int16)

        # Add the fill and scan gaps back into the results, since they may
        # have been lost
        self.logger.info('Adding fill and data gaps back into the Land'
                         ' Surface Temperature results')

        # Fix the no data locations
        no_data_locations = np.where(radiance_emitted == self.no_data_value)
        lst_data[no_data_locations] = self.no_data_value

        # Memory cleanup
        del (radiance_emitted)
        del (no_data_locations)

        product_id = self.xml_filename.split('.xml')[0]
        lst_img_filename = ''.join([product_id, '_lst', '.img'])
        lst_hdr_filename = ''.join([product_id, '_lst', '.hdr'])
        lst_aux_filename = ''.join([lst_img_filename, '.aux', '.xml'])

        self.logger.info('Creating {0}'.format(lst_img_filename))
        util.Geo.generate_raster_file(envi_driver, lst_img_filename,
                                      lst_data, x_dim, y_dim, ds_transform,
                                      ds_srs.ExportToWkt(), self.no_data_value,
                                      gdal.GDT_Int16)

        self.logger.info('Updating {0}'.format(lst_hdr_filename))
        util.Geo.update_envi_header(lst_hdr_filename, self.no_data_value)

        # Memory cleanup
        del (ds_srs)
        del (ds_transform)

        # Remove the *.aux.xml file generated by GDAL
        if os.path.exists(lst_aux_filename):
            os.unlink(lst_aux_filename)

        self.logger.info('Adding {0} to {1}'.format(lst_img_filename,
                                                    self.xml_filename))
        # Add the estimated Land Surface Temperature product to the metadata
        espa_xml = metadata_api.parse(self.xml_filename, silence=True)
        bands = espa_xml.get_bands()
        sensor_code = product_id[0:3]

        # Find the TOA Band 1 to use for the specific band details
        base_band = None
        for band in bands.band:
            if band.product == 'toa_refl' and band.name == 'toa_band1':
                base_band = band

        if base_band is None:
            raise Exception('Failed to find the TOA BLUE band'
                            ' in the input data')

        lst_band = metadata_api.band(product='lst',
                                     source='toa_refl',
                                     name='land_surface_temperature',
                                     category='image',
                                     data_type='INT16',
                                     scale_factor=SCALE_FACTOR,
                                     add_offset=0,
                                     nlines=base_band.get_nlines(),
                                     nsamps=base_band.get_nsamps(),
                                     fill_value=str(self.no_data_value))

        lst_band.set_short_name('{0}LST'.format(sensor_code))
        lst_band.set_long_name('Land Surface Temperature')
        lst_band.set_file_name(lst_img_filename)
        lst_band.set_data_units('temperature (kelvin)')

        pixel_size = metadata_api.pixel_size(base_band.pixel_size.x,
                                             base_band.pixel_size.x,
                                             base_band.pixel_size.units)
        lst_band.set_pixel_size(pixel_size)

        valid_range = metadata_api.valid_range(min=1500, max=3730)
        lst_band.set_valid_range(valid_range)

        # Set the date, but first clean the microseconds off of it
        production_date = (
            datetime.datetime.strptime(datetime.datetime.now().
                                       strftime('%Y-%m-%dT%H:%M:%S'),
                                       '%Y-%m-%dT%H:%M:%S'))

        lst_band.set_production_date(production_date)

        lst_band.set_app_version(util.Version.app_version())

        bands.add_band(lst_band)

        # Write the XML metadata file out
        with open(self.xml_filename, 'w') as fd:
            metadata_api.export(fd, espa_xml)

        # Memory cleanup
        del (lst_band)
        del (bands)
        del (espa_xml)
        del (lst_data)


# ============================================================================
if __name__ == '__main__':
    '''
    Description:
        Generate Landsat EMIS and ASTER NDVI from ASTER GED tiles for the
        specified Landsat scene.
    '''

    # Build the command line argument parser
    description = ('Reads intermediate data generated previously and combines'
                   ' them into the Land Surface Temperature product')
    parser = ArgumentParser(description=description)

    # ---- Add parameters ----
    # Required parameters
    parser.add_argument('--xml',
                        action='store', dest='xml_filename',
                        required=False, default=None,
                        help='The XML metadata file to use')

    # Optional parameters
    parser.add_argument('--version',
                        action='store_true', dest='version',
                        required=False, default=False,
                        help='Reports the version of the software')

    # Parse the command line arguments
    args = parser.parse_args()

    # Command line arguments are required so print the help if none were
    # provided
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)  # EXIT FAILURE

    # Report the version and exit
    if args.version:
        print(util.Version.version_text())
        sys.exit(0)  # EXIT SUCCESS

    # Verify that the --xml parameter was specified
    if args.xml_filename is None:
        raise Exception('--xml must be specified on the command line')
        sys.exit(1)  # EXIT FAILURE

    # Configure logging
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:%(funcName)s'
                                ' -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)

    logger = logging.getLogger(__name__)

    try:
        build_lst_data = BuildLSTData(args.xml_filename)

        # Call the main processing routine
        build_lst_data.generate_data()
    except Exception:
        logger.exception('Processing failed')
        sys.exit(1)  # EXIT FAILURE

    sys.exit(0)  # EXIT SUCCESS
