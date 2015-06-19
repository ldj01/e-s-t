#! /usr/bin/env python

'''
    FILE: l5-7_landsat_emissivity_from_aster_ged.py

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

    Date              Programmer               Reason
    ----------------  ------------------------ -------------------------------
    January/2015      Ron Dilley               Initial implementation
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
  Example HTTP file request for ASTER GED data.

      http://e4ftl01.cr.usgs.gov/ASTT/AG100.003/2000.01.01/<tile>

  Where <tile> is of the form:
      AG100.v003.44.-077.0001.h5  (Negative Longitude)
      or
      AG100.v003.44.077.0001.h5   (Positive Longitude)
'''
# Variables to hold the server name and path retrieved from the environment
SERVER_NAME = ''
SERVER_PATH = ''
# Setup some formats to apply to the URL for retrieving an ASTER GED tile
FILE_N_FORMAT = 'AG100.v003.{0:02}.{1:04}.0001'
FILE_P_FORMAT = 'AG100.v003.{0:02}.{1:03}.0001'


# ============================================================================
class EstimateLandsatEmissivity(object):
    '''
    Description:
        Defines the processor for generating the estimated Landsat emissivity
        product.
    '''

    # ------------------------------------------------------------------------
    class LandsatInfo(object):
        '''
        Description:
            Used as a data structure to cleanup parameter passing to routines.
        '''

        def __init__(self):
            self.north = None
            self.south = None
            self.east = None
            self.west = None
            self.x_pixel_size = None
            self.y_pixel_size = None
            self.min_x_extent = None
            self.max_x_extent = None
            self.min_y_extent = None
            self.max_y_extent = None
            self.dest_proj4 = None

    def __init__(self, args):
        super(EstimateLandsatEmissivity, self).__init__()

        self.xml_filename = args.xml_filename
        self.keep_intermediate_data = args.keep_intermediate_data

        # Specify the no data value we will be using, it also matches the
        # no_data_value for the ASTER data we extract and use
        self.no_data_value = -9999

    # ------------------------------------------------------------------------
    def update_envi_header(self, hdr_file_path, no_data_value):
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

    # ------------------------------------------------------------------------
    def extract_aster_data(self, emis_ds_name, ndvi_ds_name,
                           lat_ds_name, lon_ds_name):

        '''
        Description:
            Opens the requested files and extracts the associated data, which
            is then returned to the caller.

        Returns: An array of the following information.
            aster_b13_data - Aster Emissivity band 13 data.
            aster_b14_data - Aster Emissivity band 14 data.
            aster_ndvi_data - Aster NDVI band data.
            aster_lat_data - Geolocation Lattitude data.
            aster_lon_data - Geolocation Longitude data.
            x_dim - Number of columns of the data.
            y_dim - Number of rows of the data.
        '''

        try:
            # Open the Emissivity sub-dataset
            emis_sds = gdal.Open(emis_ds_name)
            if emis_sds is None:
                raise RuntimeError("GDAL failed to open {0}"
                                   .format(emis_ds_name))

            # Open the NDVI sub-dataset
            ndvi_sds = gdal.Open(ndvi_ds_name)
            if ndvi_sds is None:
                raise RuntimeError("GDAL failed to open {0}"
                                   .format(ndvi_sds_name))

            # Open the Latitude sub-dataset
            lat_sds = gdal.Open(lat_ds_name)
            if lat_sds is None:
                raise RuntimeError("GDAL failed to open {0}"
                                   .format(lat_ds_name))

            # Open the Longitude sub-dataset
            lon_sds = gdal.Open(lon_ds_name)
            if lon_sds is None:
                raise RuntimeError("GDAL failed to open {0}"
                                   .format(lon_ds_name))

            # The dimensions are the same for all the bands so just use
            # the values from the Emissivity dataset
            x_dim = emis_sds.RasterXSize
            y_dim = emis_sds.RasterYSize

            # Retrieve the band 13 data from the HDF5 input
            aster_b13_data = emis_sds.GetRasterBand(4).ReadAsArray(0, 0,
                                                                   x_dim,
                                                                   y_dim)

            # Retrieve the band 14 data from the HDF5 input
            aster_b14_data = emis_sds.GetRasterBand(5).ReadAsArray(0, 0,
                                                                   x_dim,
                                                                   y_dim)
            del (emis_sds)

            # Retrieve the NDVI band data from the HDF5 input
            aster_ndvi_data = ndvi_sds.GetRasterBand(1).ReadAsArray(0, 0,
                                                                    x_dim,
                                                                    y_dim)
            del (ndvi_sds)

            # Retrieve the Latitude data from the HDF5 input
            aster_lat_data = lat_sds.GetRasterBand(1).ReadAsArray(0, 0,
                                                                  x_dim,
                                                                  y_dim)
            del (lat_sds)

            # Retrieve the Longitude data from the HDF5 input
            aster_lon_data = lon_sds.GetRasterBand(1).ReadAsArray(0, 0,
                                                                  x_dim,
                                                                  y_dim)
            del (lon_sds)

        except Exception:
            raise

        return (aster_b13_data, aster_b14_data, aster_ndvi_data,
                aster_lat_data, aster_lon_data, x_dim, y_dim)

    # ------------------------------------------------------------------------
    def generate_raster_file(self, driver, filename, x_dim, y_dim, data,
                             geo_transform, proj_wkt, no_data_value,
                             data_type):
        '''
        Description:
            Creates a raster file on disk for the specified data, using the
            specified driver.

        Note: It is assumed that the driver supports setting of the no data
              value.
              It is the callers responsibility to fix it if it does not.

        Note: It is assumed that the caller specified the correct file
              extension in the filename parameter for the specfied driver.
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

    # ------------------------------------------------------------------------
    def mosaic_all_tiles_into_one_raster(self, src_names, dest_name):
        '''
        Description:
            Executes gdalwarp on the supplied source names to generate a
            mosaic'ed destination named file.
        '''

        logger = logging.getLogger(__name__)

        cmd = ['gdalwarp', '-wm', '2048', '-multi',
               '-srcnodata', str(self.no_data_value),
               '-dstnodata', str(self.no_data_value)]
        cmd.extend(src_names)
        cmd.append(dest_name)

        cmd = ' '.join(cmd)

        output = ''
        try:
            logger.info("Executing [{0}]".format(cmd))
            output = util.System.execute_cmd(cmd)
        except Exception:
            logger.error("Failed to mosaic tiles")
            raise
        finally:
            if len(output) > 0:
                logger.info(output)

    # ------------------------------------------------------------------------
    def warp_raster_to_match_ls_data(self, src_name, dest_name, src_proj4):
        '''
        Description:
            Executes gdalwarp on the supplied source name to generate a warped
            destination named file that matches the projection and image
            extents and pixel size specified, which were derived from Landsat
            data.
        '''

        logger = logging.getLogger(__name__)

        cmd = ['gdalwarp', '-wm', '2048', '-multi',
               '-tr', str(self.ls_info.x_pixel_size),
               str(self.ls_info.y_pixel_size),
               '-s_srs', "'" + src_proj4 + "'",
               '-t_srs', "'" + self.ls_info.dest_proj4 + "'",
               '-of', 'GTiff',
               '-overwrite', '-te',
               str(self.ls_info.min_x_extent), str(self.ls_info.min_y_extent),
               str(self.ls_info.max_x_extent), str(self.ls_info.max_y_extent),
               '-srcnodata', str(self.no_data_value),
               '-dstnodata', str(self.no_data_value)]
        cmd.append(src_name)
        cmd.append(dest_name)

        cmd = ' '.join(cmd)

        output = ''
        try:
            logger.info("Executing [{0}]".format(cmd))
            output = util.System.execute_cmd(cmd)
        except Exception:
            logger.error("Failed during warping to match Landsat")
            raise
        finally:
            if len(output) > 0:
                logger.info(output)

    # ------------------------------------------------------------------------
    def build_ls_emis_data(self, driver):
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

        # The ASTER data is in geographic projection so specify that here
        ds_srs = osr.SpatialReference()
        ds_srs.ImportFromEPSG(4326)

        # Specify the base URL to use for retrieving the ASTER GED data
        base_url_path = ''.join(['http://', SERVER_NAME, SERVER_PATH])

        # Process through the lattitude and longitude ASTER tiles which cover
        # the Landsat scene we are processing
        # - Download them
        # - Extract the Emissivity bands 13 and 14 as well as the NDVI
        # - Generate the Landsat EMIS from the 13 and 14 band data
        ls_emis_mean_filenames = list()
        aster_ndvi_mean_filenames = list()
        for lat in xrange(int(self.ls_info.south), int(self.ls_info.north)+1):
            for lon in xrange(int(self.ls_info.west),
                              int(self.ls_info.east)+1):
                # Build the base filename using the correct format
                filename = ''
                if lon < 0:
                    filename = FILE_N_FORMAT.format(lat, lon)
                else:
                    filename = FILE_P_FORMAT.format(lat, lon)

                # Build the HDF5 filename for the tile
                h5_file_path = ''.join([filename, '.h5'])

                # Build the URL and download the tile
                url_path = ''.join([base_url_path, h5_file_path])
                status_code = util.Web.http_transfer_file(url_path,
                                                          h5_file_path)

                # Check for and handle tiles that are not available in the
                # ASTER data
                if status_code != requests.codes['ok']:
                    if status_code != requests.codes['not_found']:
                        raise Exception("HTTP - Transfer Failed")
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
                    (aster_b13_data, aster_b14_data, aster_ndvi_data,
                     aster_lat_data, aster_lon_data, x_dim, y_dim) = \
                        self.extract_aster_data(emis_ds_name, ndvi_ds_name,
                                                lat_ds_name, lon_ds_name)
                except Exception:
                    logger.exception("Extracting ASTER data from tile")
                    raise

                # Remove the HDF5 tile since we have extracted all the info we
                # need from it
                if not self.keep_intermediate_data:
                    if os.path.exists(h5_file_path):
                        os.unlink(h5_file_path)

                # ------------------------------------------------------------
                # Determine the minimum and maximum latitude and longitude
                x_min = aster_lon_data.min()
                x_max = aster_lon_data.max()
                y_min = aster_lat_data.min()
                y_max = aster_lat_data.max()

                # Determine the resolution of the ASTER data
                x_res = (x_max - x_min) / float(x_dim)
                y_res = (y_max - y_min) / float(y_dim)

                # Build the geo transform to apply to the raster tile
                geo_transform = [x_min, x_res, 0, y_max, 0, -y_res]

                # ------------------------------------------------------------
                # Create a mask of the band 13 data and then scale it
                # Original matlab code is <= 0
                aster_b13_masked = np.ma.masked_less_equal(aster_b13_data, 0)
                aster_b13_scaled = aster_b13_masked * 0.001

                # Create a mask of the band 14 data and then scale it
                # Original matlab code is <= 0
                aster_b14_masked = np.ma.masked_less_equal(aster_b14_data, 0)
                aster_b14_scaled = aster_b14_masked * 0.001

                # Re-make the masks for B13 and B14 as a combined mask so the
                # estimated Landsat EMIS has the correct mask
                mask = np.ma.mask_or(np.ma.getmask(aster_b13_scaled),
                                     np.ma.getmask(aster_b14_scaled))
                aster_b13_scaled = np.ma.masked_where(mask, aster_b13_scaled)
                aster_b14_scaled = np.ma.masked_where(mask, aster_b14_scaled)

                # Memory cleanup
                del (mask)
                del (aster_b13_data)
                del (aster_b14_data)
                del (aster_b13_masked)
                del (aster_b14_masked)

                # Create the estimated Landsat EMIS data
                ls_emis_data = (0.44 * aster_b13_scaled +
                                0.4 * aster_b14_scaled + 0.156)

                # Memory cleanup
                del (aster_b13_scaled)
                del (aster_b14_scaled)

                # Create a mask of the NDVI band data and then scale it
                aster_ndvi_masked = np.ma.masked_equal(aster_ndvi_data,
                                                       self.no_data_value)
                aster_ndvi_scaled = aster_ndvi_masked * 0.01

                # Memory cleanup
                del (aster_ndvi_data)
                del (aster_ndvi_masked)

                # numpy array them so they are no longer masked
                ls_emis_numpy = np.array(ls_emis_data)
                aster_ndvi_numpy = np.array(aster_ndvi_scaled)

                # Memory cleanup
                del (ls_emis_data)
                del (aster_ndvi_scaled)

                # ------------------------------------------------------------
                # Create the estimated Landsat EMIS raster output tile
                try:
                    logger.info("Creating an estimated Landsat EMIS tile")
                    self.generate_raster_file(driver, ls_emis_tile_name,
                                              x_dim, y_dim, ls_emis_numpy,
                                              geo_transform,
                                              ds_srs.ExportToWkt(),
                                              self.no_data_value,
                                              gdal.GDT_Float32)
                except Exception:
                    logger.exception("Creating Landsat EMIS tile")
                    raise

                # ------------------------------------------------------------
                # Create the ASTER NDVI raster output tile
                try:
                    logger.info("Creating an ASTER NDVI tile")
                    self.generate_raster_file(driver, aster_ndvi_tile_name,
                                              x_dim, y_dim, aster_ndvi_numpy,
                                              geo_transform,
                                              ds_srs.ExportToWkt(),
                                              self.no_data_value,
                                              gdal.GDT_Float32)
                except Exception:
                    logger.exception("Creating ASTER NDVI tile")
                    raise

                # Memory cleanup
                del (ls_emis_numpy)
                del (aster_ndvi_numpy)

        # Check to see that we downloaded at least one ASTER tile for
        # processing.
        if len(ls_emis_mean_filenames) == 0:
            raise Exception("No ASTER tiles were downloaded")

        # Save the source proj4 string to use during warping
        src_proj4 = ds_srs.ExportToProj4()

        # Define the temporary names
        ls_emis_mosaic_name = 'landsat_emis_mosaic.tif'
        ls_emis_warped_name = 'landsat_emis_warped.tif'
        aster_ndvi_mosaic_name = 'aster_ndvi_mosaic.tif'
        aster_ndvi_warped_name = 'aster_ndvi_warped.tif'

        # Mosaic the estimated Landsat EMIS tiles into the temp EMIS
        try:
            logger.info("Building mosaic for estimated Landsat EMIS")
            self.mosaic_all_tiles_into_one_raster(ls_emis_mean_filenames,
                                                  ls_emis_mosaic_name)
        except Exception:
            logger.exception("Mosaicing EMIS tiles")
            raise

        # Mosaic the ASTER NDVI tiles into the temp NDVI
        try:
            logger.info("Building mosaic for ASTER NDVI")
            self.mosaic_all_tiles_into_one_raster(aster_ndvi_mean_filenames,
                                                  aster_ndvi_mosaic_name)
        except Exception:
            logger.exception("Mosaicing ASTER NDVI tiles")
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
            logger.info("Warping estimated Landsat EMIS to match Landsat data")
            self.warp_raster_to_match_ls_data(ls_emis_mosaic_name,
                                              ls_emis_warped_name, src_proj4)
        except Exception:
            logger.exception("Warping EMIS to match Landsat data")
            raise

        # Warp ASTER NDVI to match the Landsat data
        try:
            logger.info("Warping ASTER NDVI to match Landsat data")
            self.warp_raster_to_match_ls_data(aster_ndvi_mosaic_name,
                                              aster_ndvi_warped_name,
                                              src_proj4)
        except Exception:
            logger.exception("Warping ASTER NDVI to match Landsat data")
            raise

        if not self.keep_intermediate_data:
            # Cleanup the temp files
            if os.path.exists(ls_emis_mosaic_name):
                os.unlink(ls_emis_mosaic_name)
            if os.path.exists(aster_ndvi_mosaic_name):
                os.unlink(aster_ndvi_mosaic_name)

        return (ls_emis_warped_name, aster_ndvi_warped_name)

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

        self.ls_info = self.LandsatInfo()

        toa_bt_name = ''
        toa_green_name = ''
        toa_red_name = ''
        toa_nir_name = ''
        toa_swir1_name = ''
        toa_green_scale_factor = 1.0
        toa_red_scale_factor = 1.0
        toa_nir_scale_factor = 1.0
        toa_swir1_scale_factor = 1.0

        # Find the TOA bands to extract information from
        for band in bands.band:
            if band.product == 'toa_refl' and band.name == 'toa_band2':
                toa_green_name = band.get_file_name()
                toa_green_scale_factor = float(band.scale_factor)

            if band.product == 'toa_refl' and band.name == 'toa_band3':
                toa_red_name = band.get_file_name()
                toa_red_scale_factor = float(band.scale_factor)

            if band.product == 'toa_refl' and band.name == 'toa_band4':
                toa_nir_name = band.get_file_name()
                toa_nir_scale_factor = float(band.scale_factor)

            if band.product == 'toa_refl' and band.name == 'toa_band5':
                toa_swir1_name = band.get_file_name()
                toa_swir1_scale_factor = float(band.scale_factor)

            if band.product == 'toa_bt' and band.category == 'image':
                # Get the output pixel size
                self.ls_info.x_pixel_size = band.pixel_size.x
                self.ls_info.y_pixel_size = band.pixel_size.y

                toa_bt_name = band.get_file_name()

                # Get the output proj4 string
                self.ls_info.dest_proj4 = \
                    util.Geo.get_proj4_projection_string(toa_bt_name)

        # Error if we didn't find the required TOA bands in the data
        if len(toa_green_name) <= 0:
            raise Exception("Failed to find the TOA GREEN band"
                            " in the input data")
        if len(toa_red_name) <= 0:
            raise Exception("Failed to find the TOA RED band"
                            " in the input data")
        if len(toa_nir_name) <= 0:
            raise Exception("Failed to find the TOA NIR band"
                            " in the input data")
        if len(toa_swir1_name) <= 0:
            raise Exception("Failed to find the TOA SWIR1 band"
                            " in the input data")
        if len(toa_bt_name) <= 0:
            raise Exception("Failed to find the TOA BT band"
                            " in the input data")

        # Determine the bounding geographic coordinates for the ASTER tiles we
        # will need
        self.ls_info.north = math.ceil(gm.bounding_coordinates.north)
        self.ls_info.south = math.floor(gm.bounding_coordinates.south)
        self.ls_info.east = math.ceil(gm.bounding_coordinates.east)
        self.ls_info.west = math.floor(gm.bounding_coordinates.west)

        # Determine the UTM projection corner points
        for cp in gm.projection_information.corner_point:
            if cp.location == 'UL':
                self.ls_info.min_x_extent = cp.x
                self.ls_info.max_y_extent = cp.y
            if cp.location == 'LR':
                self.ls_info.max_x_extent = cp.x
                self.ls_info.min_y_extent = cp.y

        # Adjust the UTM coordinates for image extents becuse they are in
        # center of pixel, and we need to supply the warping with actual
        # extents
        self.ls_info.min_x_extent = (self.ls_info.min_x_extent -
                                     self.ls_info.x_pixel_size *
                                     0.5)
        self.ls_info.max_x_extent = (self.ls_info.max_x_extent +
                                     self.ls_info.x_pixel_size *
                                     0.5)
        self.ls_info.min_y_extent = (self.ls_info.min_y_extent -
                                     self.ls_info.y_pixel_size *
                                     0.5)
        self.ls_info.max_y_extent = (self.ls_info.max_y_extent +
                                     self.ls_info.y_pixel_size *
                                     0.5)

        # Save for later
        satellite = gm.satellite

        del (bands)
        del (gm)
        del (espa_xml)

        return (toa_bt_name, toa_green_name, toa_red_name,
                toa_nir_name, toa_swir1_name, toa_green_scale_factor,
                toa_red_scale_factor, toa_nir_scale_factor,
                toa_swir1_scale_factor, satellite)

    # ------------------------------------------------------------------------
    # TODO - NEED TO PROCESS A COASTAL SCENE BECAUSE MAY NOT HAVE ASTER TILE
    #        DATA FOR SOME OF THE SCENE AND WILL NEED TO DO SOMETHING ELSE
    #      - I have implemented a solution for this, but have not been able to
    #        test it.  I have not been able to find a scene that has this
    #        condition.
    # ------------------------------------------------------------------------
    def generate_product(self):
        '''
        Description:
            Provides the main processing algorithm for generating the
            estimated Landsat emissivity product.  It produces the final
            emissivity product.
        '''

        logger = logging.getLogger(__name__)

        try:
            (toa_bt_name, toa_green_name, toa_red_name, toa_nir_name,
             toa_swir1_name, toa_green_scale_factor, toa_red_scale_factor,
             toa_nir_scale_factor, toa_swir1_scale_factor,
             satellite) = self.retrieve_metadata_information()
        except:
            logger.exception("Failed reading input XML metadata file")
            raise

        # Register all the gdal drivers and choose the GeoTiff for our temp
        # output
        gdal.AllRegister()
        geotiff_driver = gdal.GetDriverByName('GTiff')
        envi_driver = gdal.GetDriverByName('ENVI')

        # Read the Landsat bands into memory
        logger.info("Loading Landsat TOA input bands")
        # GREEN
        ds = gdal.Open(toa_green_name)
        x_dim = ds.RasterXSize  # They are all the same size
        y_dim = ds.RasterYSize
        ls_green_data = ds.GetRasterBand(1).ReadAsArray(0, 0, x_dim, y_dim)
        ls_green_masked = np.ma.masked_equal(ls_green_data,
                                             self.no_data_value)
        ls_green_masked = ls_green_masked * toa_green_scale_factor

        # RED
        ds = gdal.Open(toa_red_name)
        ls_red_data = ds.GetRasterBand(1).ReadAsArray(0, 0, x_dim, y_dim)
        ls_red_masked = np.ma.masked_equal(ls_red_data, self.no_data_value)
        ls_red_masked = ls_red_masked * toa_red_scale_factor

        # NIR
        ds = gdal.Open(toa_nir_name)
        ls_nir_data = ds.GetRasterBand(1).ReadAsArray(0, 0, x_dim, y_dim)
        ls_nir_masked = np.ma.masked_equal(ls_nir_data, self.no_data_value)
        ls_nir_masked = ls_nir_masked * toa_nir_scale_factor

        # SWIR1
        ds = gdal.Open(toa_swir1_name)
        ls_swir1_data = ds.GetRasterBand(1).ReadAsArray(0, 0, x_dim, y_dim)
        ls_swir1_masked = np.ma.masked_equal(ls_swir1_data,
                                             self.no_data_value)
        ls_swir1_masked = ls_swir1_masked * toa_swir1_scale_factor

        # Save the locations of the fill and scan gaps
        ls_green_no_data_locations = \
            np.where(ls_green_data == self.no_data_value)
        ls_red_no_data_locations = \
            np.where(ls_red_data == self.no_data_value)
        ls_nir_no_data_locations = \
            np.where(ls_nir_data == self.no_data_value)
        ls_swir1_no_data_locations = \
            np.where(ls_swir1_data == self.no_data_value)

        # Save for the output products
        ds_tmp_srs = osr.SpatialReference()
        ds_tmp_srs.ImportFromWkt(ds.GetProjection())
        ds_tmp_transform = ds.GetGeoTransform()

        # Memory cleanup
        del (ls_green_data)
        del (ls_red_data)
        del (ls_nir_data)
        del (ls_swir1_data)
        del (ds)

        # Re-make the masks for NIR and RED as a combined mask so NDVI has the
        # correct mask
        mask = np.ma.mask_or(np.ma.getmask(ls_nir_masked),
                             np.ma.getmask(ls_red_masked))
        ls_nir_masked = np.ma.masked_where(mask, ls_nir_masked)
        ls_red_masked = np.ma.masked_where(mask, ls_red_masked)

        # Re-make the masks for GREEN and SWIR1 as a combined mask so NDSI has
        # the correct mask
        mask = np.ma.mask_or(np.ma.getmask(ls_green_masked),
                             np.ma.getmask(ls_swir1_masked))
        ls_green_masked = np.ma.masked_where(mask, ls_green_masked)
        ls_swir1_masked = np.ma.masked_where(mask, ls_swir1_masked)

        # Build the Landsat TOA NDVI data
        logger.info("Building TOA based NDVI for Landsat data")
        ls_ndvi_masked = ((ls_nir_masked - ls_red_masked) /
                          (ls_nir_masked + ls_red_masked))

        # Build the Landsat TOA NDSI data
        logger.info("Building TOA based NDSI for Landsat data")
        ls_ndsi_masked = ((ls_green_masked - ls_swir1_masked) /
                          (ls_green_masked + ls_swir1_masked))

        # Memory cleanup
        del (mask)
        del (ls_red_masked)
        del (ls_nir_masked)
        del (ls_green_masked)
        del (ls_swir1_masked)

        # Harden the mask so it does not change in the following statement
        ls_ndvi_masked = ls_ndvi_masked.harden_mask()
        # Turn all negative values to zero
        ls_ndvi_masked[ls_ndvi_masked < 0] = 0

        # Save the locations for the specfied snow pixels
        logger.info("Determine snow pixel locations")
        selected_snow_locations = np.where(ls_ndsi_masked > 0.4)

        # Memory cleanup
        del (ls_ndsi_masked)

        # Build the estimated Landsat EMIS data from the ASTER GED data and
        # warp it to the Landsat scenes projection and image extents
        # For convenience the ASTER NDVI is also extracted and warped to the
        # Landsat scenes projection and image extents
        logger.info("Build thermal emissivity band and retrieve ASTER NDVI")
        (ls_emis_warped_name,
         aster_ndvi_warped_name) = self.build_ls_emis_data(geotiff_driver)

        # Load the warped estimated Landsat EMIS into memory
        ds = gdal.Open(ls_emis_warped_name)
        ls_emis_data = ds.GetRasterBand(1).ReadAsArray(0, 0, x_dim, y_dim)
        ls_emis_masked = np.ma.masked_equal(ls_emis_data, self.no_data_value)
        # Harden the mask so it does not change
        ls_emis_masked = ls_emis_masked.harden_mask()

        # Load the warped ASTER NDVI into memory
        ds = gdal.Open(aster_ndvi_warped_name)
        aster_ndvi_data = ds.GetRasterBand(1).ReadAsArray(0, 0, x_dim, y_dim)
        aster_ndvi_masked = np.ma.masked_equal(aster_ndvi_data,
                                               self.no_data_value)
        # Harden the mask so it does not change
        aster_ndvi_masked = aster_ndvi_masked.harden_mask()

        # Save the locations of the gaps in the estimated Landsat EMIS
        # and Aster NDVI data
        ls_emis_gap_locations = np.where(ls_emis_data == 0)
        ls_emis_no_data_locations = \
            np.where(ls_emis_data == self.no_data_value)
        aster_ndvi_gap_locations = np.where(aster_ndvi_data == 0)
        aster_ndvi_no_data_locations = np.where(aster_ndvi_data ==
                                                self.no_data_value)

        # Memory cleanup
        del (ds)
        del (ls_emis_data)
        del (aster_ndvi_data)

        if not self.keep_intermediate_data:
            # Cleanup the temp files since we have them in memory
            if os.path.exists(ls_emis_warped_name):
                os.unlink(ls_emis_warped_name)
            if os.path.exists(aster_ndvi_warped_name):
                os.unlink(aster_ndvi_warped_name)

        logger.info("Normalizing Landsat and ASTER NDVI")
        # Normalize Landsat NDVI by max value
        min_ls_ndvi = ls_ndvi_masked.min()
        max_ls_ndvi = ls_ndvi_masked.max()
        ls_ndvi_masked = ls_ndvi_masked / float(max_ls_ndvi)

        # Normalize ASTER NDVI by max value
        min_aster_ndvi = aster_ndvi_masked.min()
        max_aster_ndvi = aster_ndvi_masked.max()
        aster_ndvi_masked = aster_ndvi_masked / float(max_aster_ndvi)

        # Calculate fractional veg-cover for both Landsat and ASTER NDVI
        logger.info("Calculating fractional vegetation cover")
        fv_Landsat = 1.0 - ((max_ls_ndvi - ls_ndvi_masked) /
                            (max_ls_ndvi - min_ls_ndvi))
        fv_Aster = 1.0 - ((max_aster_ndvi - aster_ndvi_masked) /
                          (max_aster_ndvi - min_aster_ndvi))

        # Memory cleanup
        del (ls_ndvi_masked)
        del (aster_ndvi_masked)

        ls_soil = None
        ls_mod = None
        ls_emis_final = None
        # Adjust estimated Landsat EMIS for vegetation and snow, to generate
        # the final Landsat EMIS data
        if satellite == 'LANDSAT_7':
            logger.info("Adjusting estimated Landsat 7 EMIS"
                        " for vegetation and snow")
            ls_soil = ((ls_emis_masked - 0.975 * fv_Aster) / (1.0 - fv_Aster))
            ls_mod = (0.9848 * fv_Landsat + ls_soil * (1.0 - fv_Landsat))

            # Create a copy as a non-masked array
            ls_emis_final = np.array(ls_mod)
            # Medium snow
            ls_emis_final[selected_snow_locations] = 0.9869

        elif satellite == 'LANDSAT_5':
            logger.info("Adjusting estimated Landsat 5 EMIS"
                        " for vegetation and snow")
            ls_soil = ((ls_emis_masked - 0.975 * fv_Aster) / (1.0 - fv_Aster))
            ls_mod = (0.9851 * fv_Landsat + ls_soil * (1.0 - fv_Landsat))

            # Create a copy as a non-masked array
            ls_emis_final = np.array(ls_mod)
            # Medium snow
            ls_emis_final[selected_snow_locations] = 0.9851

        # Memory cleanup
        del (fv_Landsat)
        del (fv_Aster)
        del (ls_soil)
        del (ls_mod)
        del (selected_snow_locations)

        # Add the fill and scan gaps and ASTER gaps back into the results,
        # since they may have been lost
        logger.info("Adding fill and data gaps back into the"
                    " Landsat emissivity results")
        ls_emis_final[ls_emis_no_data_locations] = self.no_data_value
        ls_emis_final[ls_emis_gap_locations] = self.no_data_value
        ls_emis_final[aster_ndvi_no_data_locations] = self.no_data_value
        ls_emis_final[aster_ndvi_gap_locations] = self.no_data_value
        ls_emis_final[ls_green_no_data_locations] = self.no_data_value
        ls_emis_final[ls_red_no_data_locations] = self.no_data_value
        ls_emis_final[ls_nir_no_data_locations] = self.no_data_value
        ls_emis_final[ls_swir1_no_data_locations] = self.no_data_value

        # Memory cleanup
        del (ls_emis_no_data_locations)
        del (ls_emis_gap_locations)
        del (aster_ndvi_no_data_locations)
        del (aster_ndvi_gap_locations)
        del (ls_green_no_data_locations)
        del (ls_red_no_data_locations)
        del (ls_nir_no_data_locations)
        del (ls_swir1_no_data_locations)

        product_id = self.xml_filename.split('.xml')[0]
        ls_emis_img_filename = ''.join([product_id, '_emis', '.img'])
        ls_emis_hdr_filename = ''.join([product_id, '_emis', '.hdr'])
        ls_emis_aux_filename = ''.join([ls_emis_img_filename, '.aux', '.xml'])

        logger.info("Creating {0}".format(ls_emis_img_filename))
        self.generate_raster_file(envi_driver, ls_emis_img_filename,
                                  x_dim, y_dim, ls_emis_final,
                                  ds_tmp_transform, ds_tmp_srs.ExportToWkt(),
                                  self.no_data_value, gdal.GDT_Float32)

        logger.info("Updating {0}".format(ls_emis_hdr_filename))
        self.update_envi_header(ls_emis_hdr_filename, self.no_data_value)

        # Remove the *.aux.xml file generated by GDAL
        if os.path.exists(ls_emis_aux_filename):
            os.unlink(ls_emis_aux_filename)

        # Add the estimated Landsat emissivity to the metadata XML
        espa_xml = metadata_api.parse(self.xml_filename, silence=True)
        bands = espa_xml.get_bands()
        sensor_code = product_id[0:3]
        source_product = 'toa_refl'

        # Find the TOA Band 1 to use for the specific band details
        base_band = None
        for band in bands.band:
            if band.product == source_product and band.name == 'toa_band1':
                base_band = band

        if base_band is None:
            raise Exception("Failed to find the TOA BLUE band"
                            " in the input data")

        emis_band = metadata_api.band(product="lst_temp",
                                      source=source_product,
                                      name="landsat_emis",
                                      category="image",
                                      data_type="FLOAT32",
                                      nlines=base_band.get_nlines(),
                                      nsamps=base_band.get_nsamps(),
                                      fill_value=str(self.no_data_value))

        emis_band.set_short_name('{0}EMIS'.format(sensor_code))
        emis_band.set_long_name("Landsat emissivity estimated from ASTER GED"
                                " data")
        emis_band.set_file_name(ls_emis_img_filename)
        emis_band.set_data_units("Emissivity Coefficient")

        pixel_size = metadata_api.pixel_size(base_band.pixel_size.x,
                                             base_band.pixel_size.x,
                                             base_band.pixel_size.units)
        emis_band.set_pixel_size(pixel_size)

        valid_range = metadata_api.valid_range(min=0.0, max=1.0)
        emis_band.set_valid_range(valid_range)

        # Set the date, but first clean the microseconds off of it
        production_date = \
            datetime.datetime.strptime(datetime.datetime.now().
                                       strftime('%Y-%m-%dT%H:%M:%S'),
                                       '%Y-%m-%dT%H:%M:%S')

        emis_band.set_production_date(production_date)

        emis_band.set_app_version(util.Version.app_version())

        bands.add_band(emis_band)

        # Write the XML metadata file out
        with open(self.xml_filename, 'w') as fd:
            metadata_api.export(fd, espa_xml)

        # Memory cleanup
        del (ls_emis_final)


# ============================================================================
if __name__ == '__main__':
    '''
    Description:
        Generate Landsat EMIS and ASTER NDVI from ASTER GED tiles for the
        specified Landsat scene.
    '''

    # Build the command line argument parser
    description = ("Retrieves ASTER GED data and derives an estimated Landsat"
                   "emissivity product for the specified Landsat scene")
    parser = ArgumentParser(description=description)

    # ---- Add parameters ----
    # Required parameters
    parser.add_argument('--xml',
                        action='store', dest='xml_filename',
                        required=False, default=None,
                        help="The XML metadata file to use")

    # Optional parameters
    parser.add_argument('--keep-intermediate-data',
                        action='store_true', dest='keep_intermediate_data',
                        required=False, default=False,
                        help="Keep any intermediate data generated")

    parser.add_argument('--debug',
                        action='store_true', dest='debug',
                        required=False, default=False,
                        help="Turn debug messaging on")

    parser.add_argument('--version',
                        action='store_true', dest='version',
                        required=False, default=False,
                        help="Reports the version of the software")

    # Parse the command line arguments
    args = parser.parse_args()

    # Report the version and exit
    if args.version:
        print "Version: {0}".format(util.Version.version_text())
        sys.exit(0)  # EXIT SUCCESS

    # Verify that the --xml parameter was specified
    if args.xml_filename is None:
        raise Exception("--xml must be specified on the command line")
        sys.exit(1)  # EXIT FAILURE

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

    try:
        # Grab the server name from the environment
        if 'ASTER_GED_SERVER_NAME' not in os.environ:
            raise Exception("Environment variable ASTER_GED_SERVER_NAME is"
                            " not defined")
        else:
            SERVER_NAME = os.environ.get('ASTER_GED_SERVER_NAME')

        # Grab the server path from the environment or default it
        SERVER_PATH = os.environ.get('ASTER_GED_SERVER_PATH',
                                     '/ASTT/AG100.003/2000.01.01/')

        processor = EstimateLandsatEmissivity(args)

        # Call the main processing routine
        processor.generate_product()
    except Exception, e:
        if hasattr(e, 'output'):
            logger.error("Output [%s]" % e.output)
        logger.exception("Processing failed")
        sys.exit(1)  # EXIT FAILURE

    sys.exit(0)  # EXIT SUCCESS
