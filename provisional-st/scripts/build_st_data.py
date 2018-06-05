#! /usr/bin/env python

'''
    FILE: build_st_data.py

    PURPOSE: Builds the ST product from the intermediate data that was
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
'''

import os
import sys
import logging
import datetime
from argparse import ArgumentParser
from lxml import objectify as objectify
from osgeo import gdal, osr
import numpy as np
from espa import Metadata

# Import local modules
import st_utilities as util


SCALE_FACTOR = 0.1
MULT_FACTOR = 10.0


class BuildSTData(object):
    '''
    Description:
        Defines the processor for generating the Surface Temperature product.
    '''

    def __init__(self, xml_filename):
        super(BuildSTData, self).__init__()

        # Keep local copies of this
        self.xml_filename = xml_filename

        self.st_data_dir = ''
        # Grab the data directory from the environment
        if 'ST_DATA_DIR' not in os.environ:
            raise Exception('Environment variable ST_DATA_DIR is'
                            ' not defined')
        else:
            self.st_data_dir = os.environ.get('ST_DATA_DIR')

        # Specify the no data value we will be using, it also matches the
        # no_data_value for the ASTER data we extract and use
        self.no_data_value = -9999

        # Setup the logger to use
        self.logger = logging.getLogger(__name__)

        # Setup names
        self.thermal_name = ''
        self.transmittance_name = ''
        self.upwelled_name = ''
        self.downwelled_name = ''
        self.emissivity_name = ''
        self.satellite = ''

    def retrieve_metadata_information(self):
        '''
        Description:
            Loads and reads required information from the metadata XML file.
        '''

        # Read the XML metadata
        metadata = Metadata(xml_filename=self.xml_filename)

        self.thermal_name = ''
        self.transmittance_name = ''
        self.upwelled_name = ''
        self.downwelled_name = ''
        self.emissivity_name = ''

        # Find the TOA bands to extract information from
        for band in metadata.xml_object.bands.band:
            if (band.get("product") == 'st_intermediate' and
                    band.get("name") == 'st_thermal_radiance'):
                self.thermal_name = str(band.file_name)

            if (band.get("product") == 'st_intermediate' and
                    band.get("name") == 'st_atmospheric_transmittance'):
                self.transmittance_name = str(band.file_name)

            if (band.get("product") == 'st_intermediate' and
                    band.get("name") == 'st_upwelled_radiance'):
                self.upwelled_name = str(band.file_name)

            if (band.get("product") == 'st_intermediate' and
                    band.get("name") == 'st_downwelled_radiance'):
                self.downwelled_name = str(band.file_name)

            if (band.get("product") == 'st_intermediate' and
                    band.get("name") == 'emis'):
                self.emissivity_name = str(band.file_name)

        # Error if we didn't find the required TOA bands in the data
        if len(self.thermal_name) <= 0:
            raise Exception('Failed to find the st_thermal_radiance band'
                            ' in the input data')
        if len(self.transmittance_name) <= 0:
            raise Exception('Failed to find the st_atmospheric_transmittance'
                            ' in the input data')
        if len(self.upwelled_name) <= 0:
            raise Exception('Failed to find the st_upwelled_radiance'
                            ' in the input data')
        if len(self.downwelled_name) <= 0:
            raise Exception('Failed to find the st_downwelled_radiance'
                            ' in the input data')
        if len(self.emissivity_name) <= 0:
            raise Exception('Failed to find the emis'
                            ' in the input data')

        # Save for later
        self.satellite = metadata.xml_object.global_metadata.satellite

        del metadata

    def generate_data(self):
        '''
        Description:
            Provides the main processing algorithm for building the Surface 
            Temperature product.  It produces the final ST product.
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
        self.logger.info('Loading intermediate thermal band data [{0}]'
                         .format(self.thermal_name))
        dataset = gdal.Open(self.thermal_name)
        x_dim = dataset.RasterXSize  # They are all the same size
        y_dim = dataset.RasterYSize

        thermal_data = dataset.GetRasterBand(1).ReadAsArray(0, 0, x_dim, y_dim)

        # Atmospheric transmittance
        self.logger.info('Loading intermediate transmittance band data [{0}]'
                         .format(self.transmittance_name))
        dataset = gdal.Open(self.transmittance_name)
        trans_data = dataset.GetRasterBand(1).ReadAsArray(0, 0, x_dim, y_dim)

        # Atmospheric path radiance - upwelled radiance
        self.logger.info('Loading intermediate upwelled band data [{0}]'
                         .format(self.upwelled_name))
        dataset = gdal.Open(self.upwelled_name)
        upwelled_data = dataset.GetRasterBand(1).ReadAsArray(0, 0, x_dim, y_dim)

        self.logger.info('Calculating surface radiance')
        # Surface radiance
        with np.errstate(invalid='ignore'):
            surface_radiance = (thermal_data - upwelled_data) / trans_data

        # Fix the no data locations
        no_data_locations = np.where(thermal_data == self.no_data_value)
        surface_radiance[no_data_locations] = self.no_data_value

        no_data_locations = np.where(trans_data == self.no_data_value)
        surface_radiance[no_data_locations] = self.no_data_value

        no_data_locations = np.where(upwelled_data == self.no_data_value)
        surface_radiance[no_data_locations] = self.no_data_value

        # Memory cleanup
        del thermal_data
        del trans_data
        del upwelled_data
        del no_data_locations

        # Downwelling sky irradiance
        self.logger.info('Loading intermediate downwelled band data [{0}]'
                         .format(self.downwelled_name))
        dataset = gdal.Open(self.downwelled_name)
        downwelled_data = dataset.GetRasterBand(1).ReadAsArray(0, 0, x_dim,
                                                               y_dim)

        # Landsat emissivity estimated from ASTER GED data
        self.logger.info('Loading intermediate emissivity band data [{0}]'
                         .format(self.emissivity_name))
        dataset = gdal.Open(self.emissivity_name)
        emissivity_data = dataset.GetRasterBand(1).ReadAsArray(0, 0, x_dim,
                                                               y_dim)

        # Save for the output product
        ds_srs = osr.SpatialReference()
        ds_srs.ImportFromWkt(dataset.GetProjection())
        ds_transform = dataset.GetGeoTransform()

        # Memory cleanup
        del dataset

        # Estimate Earth-emitted radiance by subtracting off the reflected
        # downwelling component
        radiance = (surface_radiance -
                    (1.0 - emissivity_data) * downwelled_data)

        # Account for surface emissivity to get Plank emitted radiance
        self.logger.info('Calculating Plank emitted radiance')
        with np.errstate(invalid='ignore'):
            radiance_emitted = radiance / emissivity_data

        # Fix the no data locations
        no_data_locations = np.where(surface_radiance == self.no_data_value)
        radiance_emitted[no_data_locations] = self.no_data_value

        no_data_locations = np.where(downwelled_data == self.no_data_value)
        radiance_emitted[no_data_locations] = self.no_data_value

        no_data_locations = np.where(emissivity_data == self.no_data_value)
        radiance_emitted[no_data_locations] = self.no_data_value

        # Memory cleanup
        del downwelled_data
        del emissivity_data
        del surface_radiance
        del radiance
        del no_data_locations

        # Use Brightness Temperature LUT to get skin temperature
        # Read the correct one for what we are processing
        if self.satellite == 'LANDSAT_8':
            self.logger.info('Using Landsat 8 Brightness Temperature LUT')
            bt_name = 'L8_Brightness_Temperature_LUT.txt'

        elif self.satellite == 'LANDSAT_7':
            self.logger.info('Using Landsat 7 Brightness Temperature LUT')
            bt_name = 'L7_Brightness_Temperature_LUT.txt'

        elif self.satellite == 'LANDSAT_5':
            self.logger.info('Using Landsat 5 Brightness Temperature LUT')
            bt_name = 'L5_Brightness_Temperature_LUT.txt'

        elif self.satellite == 'LANDSAT_4':
            self.logger.info('Using Landsat 4 Brightness Temperature LUT')
            bt_name = 'L4_Brightness_Temperature_LUT.txt'

        bt_data = np.loadtxt(os.path.join(self.st_data_dir, bt_name),
                             dtype=float, delimiter=' ')
        bt_radiance_lut = bt_data[:, 1]
        bt_temp_lut = bt_data[:, 0]

        self.logger.info('Generating ST results')
        st_data = np.interp(radiance_emitted, bt_radiance_lut, bt_temp_lut)

        # Scale the result
        st_data = st_data * MULT_FACTOR

        # Add the fill and scan gaps back into the results, since they may
        # have been lost
        self.logger.info('Adding fill and data gaps back into the Surface'
                         ' Temperature results')

        # Fix the no data locations
        no_data_locations = np.where(radiance_emitted == self.no_data_value)
        st_data[no_data_locations] = self.no_data_value

        # Memory cleanup
        del radiance_emitted
        del no_data_locations

        product_id = self.xml_filename.split('.xml')[0]
        st_img_filename = ''.join([product_id, '_st', '.img'])
        st_hdr_filename = ''.join([product_id, '_st', '.hdr'])
        st_aux_filename = ''.join([st_img_filename, '.aux', '.xml'])

        self.logger.info('Creating {0}'.format(st_img_filename))
        util.Geo.generate_raster_file(envi_driver, st_img_filename,
                                      st_data, x_dim, y_dim, ds_transform,
                                      ds_srs.ExportToWkt(), self.no_data_value,
                                      gdal.GDT_Int16)

        self.logger.info('Updating {0}'.format(st_hdr_filename))
        util.Geo.update_envi_header(st_hdr_filename, self.no_data_value)

        # Memory cleanup
        del ds_srs
        del ds_transform

        # Remove the *.aux.xml file generated by GDAL
        if os.path.exists(st_aux_filename):
            os.unlink(st_aux_filename)

        self.logger.info('Adding {0} to {1}'.format(st_img_filename,
                                                    self.xml_filename))
        # Add the estimated Surface Temperature product to the metadata
        metadata = Metadata(xml_filename=self.xml_filename)

        # Create an element maker and band element
        em = objectify.ElementMaker(annotate=False, namespace=None, nsmap=None)
        st_band = em.band()

        # Find the TOA Band 1 to use for the specific band details
        base_band = None
        for band in metadata.xml_object.bands.band:
            if (band.get("product") == 'toa_refl' and
                    band.get("name") == 'toa_band1'):
                base_band = band

        if base_band is None:
            raise Exception('Failed to find the TOA band 1'
                            ' in the input data')

        # Set attributes for the band element
        st_band.set('product', 'st')
        st_band.set('source', 'toa_refl')
        st_band.set('name', 'surface_temperature')
        st_band.set('category', 'image')
        st_band.set('data_type', 'INT16')
        st_band.set('scale_factor', "%.6f" % SCALE_FACTOR)
        st_band.set('add_offset', "%.6f" % 0.0)
        st_band.set('nlines', str(int(base_band.get("nlines"))))
        st_band.set('nsamps', str(int(base_band.get("nsamps"))))
        st_band.set('fill_value', str(self.no_data_value))

        # Add elements to the band object
        st_band.short_name = em.element('{0}ST'.format(product_id[0:4]))
        st_band.long_name = em.element('Surface Temperature')
        st_band.file_name = em.element(str(st_img_filename))

        # Create a pixel size element and add attributes to it
        st_band.pixel_size = em.element()
        st_band.pixel_size.set('x', str(base_band.pixel_size.get('x')))
        st_band.pixel_size.set('y', str(base_band.pixel_size.get('x')))
        st_band.pixel_size.set('units', str(base_band.pixel_size.get('units')))

        st_band.resample_method = em.element('none')
        st_band.data_units = em.element('temperature (kelvin)')

        # Create a valid range element and add attributes to it
        st_band.valid_range = em.element()
        st_band.valid_range.set('min', '1500')
        st_band.valid_range.set('max', '3730')

        st_band.app_version = em.element(str(util.Version.app_version()))

        # Set the date, but first clean the microseconds off of it
        production_date = ('{0}Z'.format(datetime.datetime.now()
                                         .strftime('%Y-%m-%dT%H:%M:%S')))
        st_band.production_date = em.element(str(production_date))

        # Add the new band to the XML, validate it, and write it
        metadata.xml_object.bands.append(st_band)
        metadata.validate()
        metadata.write(xml_filename=self.xml_filename)

        # Memory cleanup
        del metadata
        del st_band


def main():
    '''
    Description:
        Generate Landsat EMIS and ASTER NDVI from ASTER GED tiles for the
        specified Landsat scene.
    '''

    # Build the command line argument parser
    parser = ArgumentParser(description='Reads intermediate data generated'
                                        ' previously and combines them into'
                                        ' the Surface Temperature product')

    # ---- Add parameters ----
    # Required parameters
    parser.add_argument('--xml',
                        action='store', dest='xml_filename',
                        required=False, default=None,
                        help='The XML metadata file to use')

    # Optional parameters
    parser.add_argument('--version',
                        action='version',
                        version=util.Version.version_text())

    # Parse the command line arguments
    args = parser.parse_args()

    # Command line arguments are required so print the help if none were
    # provided
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)  # EXIT FAILURE

    # Verify that the --xml parameter was specified
    if args.xml_filename is None:
        raise Exception('--xml must be specified on the command line')

    # Configure logging
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:%(funcName)s'
                                ' -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.DEBUG)

    logger = logging.getLogger(__name__)

    try:
        build_st_data = BuildSTData(args.xml_filename)

        # Call the main processing routine
        build_st_data.generate_data()
    except Exception:
        logger.exception('Processing failed')
        sys.exit(1)  # EXIT FAILURE

    sys.exit(0)  # EXIT SUCCESS

if __name__ == '__main__':
    main()
