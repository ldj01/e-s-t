#! /usr/bin/env python

import os
import sys
import glob
import logging
import math
import requests
import commands
import numpy as np
from cStringIO import StringIO
from argparse import ArgumentParser
from osgeo import gdal, osr
from time import sleep

# espa-common objects and methods
from espa_constants import EXIT_FAILURE
from espa_constants import EXIT_SUCCESS
import metadata_api


'''
  Example HTTP file request for ASTER GED data.
      http://e4ftl01.cr.usgs.gov/ASTT/AG100.003/2000.01.01/AG100.v003.44.-077.0001.h5
'''
# Setup some formats to generate the URL to retrieve an ASTER GED tile
SERVER = 'http://e4ftl01.cr.usgs.gov'
SERVER_PATH = '/ASTT/AG100.003/2000.01.01/'
FILE_N_FORMAT = 'AG100.v003.{0:02}.{1:04}.0001'
FILE_P_FORMAT = 'AG100.v003.{0:02}.{1:03}.0001'


# Specify the no data value we will be using, it also matches the
# no_data_value for the ASTER data we extract and use
NO_DATA_VALUE = -9999


# ============================================================================
def execute_cmd(cmd):
    '''
    Description:
      Execute a command line and return the terminal output or raise an
      exception

    Returns:
        output - The stdout and/or stderr from the executed command.
    '''

    output = ''

    (status, output) = commands.getstatusoutput(cmd)

    if status < 0:
        message = "Application terminated by signal [%s]" % cmd
        if len(output) > 0:
            message = ' Stdout/Stderr is: '.join([message, output])
        raise Exception(message)

    if status != 0:
        message = "Application failed to execute [%s]" % cmd
        if len(output) > 0:
            message = ' Stdout/Stderr is: '.join([message, output])
        raise Exception(message)

    if os.WEXITSTATUS(status) != 0:
        message = "Application [%s] returned error code [%d]" \
                  % (cmd, os.WEXITSTATUS(status))
        if len(output) > 0:
            message = ' Stdout/Stderr is: '.join([message, output])
        raise Exception(message)

    return output


# ============================================================================
def http_transfer_file(download_url, destination_file):
    '''
    Description:
      Using http transfer a file from a source location to a destination
      file on the localhost.
    '''

    logger = logging.getLogger(__name__)

    logger.info(download_url)

    session = requests.Session()

    session.mount('http://', requests.adapters.HTTPAdapter(max_retries=3))
    session.mount('https://', requests.adapters.HTTPAdapter(max_retries=3))

    retry_attempt = 0
    done = False
    while not done:
        req = None
        try:
            req = session.get(url=download_url, timeout=300.0)

            if not req.ok:
                logger.error("Transfer Failed - HTTP")
                req.raise_for_status()

            with open(destination_file, 'wb') as local_fd:
                local_fd.write(req.content)

            done = True

        except:
            logger.exception("Transfer Issue - HTTP")
            if retry_attempt > 3:
                raise Exception("Transfer Failed - HTTP"
                                " - exceeded retry limit")
            retry_attempt += 1
            sleep(int(1.5 * retry_attempt))

        finally:
            if req is not None:
                req.close()

    logger.info("Transfer Complete - HTTP")


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
            if (line.startswith('data ignore value')
                    or line.startswith('description')):
                pass
            else:
                sb.write(line)

            if line.startswith('description'):
                # This may be on multiple lines so read lines until
                # we find the closing brace
                if not line.strip().endswith('}'):
                    while 1:
                        next_line = tmp_fd.readline()
                        if (not next_line
                                or next_line.strip().endswith('}')):
                            break
                sb.write('description = {ESPA-generated file}\n')
            elif (line.startswith('data type')
                  and (no_data_value is not None)):
                sb.write('data ignore value = %s\n' % no_data_value)

    # Do the actual replace here
    with open(hdr_file_path, 'w') as tmp_fd:
        tmp_fd.write(sb.getvalue())


# ============================================================================
def extract_aster_data(emis_ds_name, lat_ds_name, lon_ds_name):

    '''
    Description:
        Opens the requested files and extracts the associated data, which is
        then returned to the caller.

    Returns: An array of the following information.
        Aster Emissivity band 13 data.
        Aster Emissivity band 14 data.
        Aster NDVI band data.
        Geolocation Lattitude data.
        Geolocation Longitude data.
        Number of columns (x_dim) of the data.
        Number of rows (y_dim) of the data.
    '''

    try:
        # Open the Emissivity sub-dataset
        emis_sds = gdal.Open(emis_ds_name)
        if emis_sds is None:
            raise RuntimeError("GDAL failed to open {0}".format(emis_ds_name))

        # Open the Latitude sub-dataset
        lat_sds = gdal.Open(lat_ds_name)
        if lat_sds is None:
            raise RuntimeError("GDAL failed to open {0}".format(lat_ds_name))

        # Open the Longitude sub-dataset
        lon_sds = gdal.Open(lon_ds_name)
        if lon_sds is None:
            raise RuntimeError("GDAL failed to open {0}".format(lon_ds_name))

        # The dimensions are the same for all the bands so just use
        # the values from the Emissivity dataset
        x_dim = emis_sds.RasterXSize
        y_dim = emis_sds.RasterYSize

        # Retrieve the band 13 data from the HDF5 input
        ds_13_band_data = emis_sds.GetRasterBand(4).ReadAsArray(0, 0,
                                                                x_dim, y_dim)

        # Retrieve the band 14 data from the HDF5 input
        ds_14_band_data = emis_sds.GetRasterBand(5).ReadAsArray(0, 0,
                                                                x_dim, y_dim)
        del (emis_sds)

        # Retrieve the Latitude data from the HDF5 input
        ds_lat_band_data = lat_sds.GetRasterBand(1).ReadAsArray(0, 0,
                                                                x_dim, y_dim)
        del (lat_sds)

        # Retrieve the Longitude data from the HDF5 input
        ds_lon_band_data = lon_sds.GetRasterBand(1).ReadAsArray(0, 0,
                                                                x_dim, y_dim)
        del (lon_sds)

    except Exception:
        raise

    return (ds_13_band_data, ds_14_band_data,
            ds_lat_band_data, ds_lon_band_data, x_dim, y_dim)


# ============================================================================
def generate_raster_file(driver, filename, x_dim, y_dim, data,
                         geo_transform, proj_wkt, no_data_value, data_type):
    '''
    Description:
        Creates a raster file on disk for the specified data, using the
        specified driver.

    Note: It is assumed that the driver supports setting of the no data value.
          It is the callers responsibility to fix it if it does not.
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
def get_original_projection(img_filename):

    ds = gdal.Open(img_filename)
    if ds is None:
        raise RuntimeError("GDAL failed to open (%s)" % img_filename)

    ds_srs = osr.SpatialReference()
    ds_srs.ImportFromWkt(ds.GetProjection())

    proj4 = ds_srs.ExportToProj4()

    del (ds_srs)
    del (ds)

    return proj4


# ============================================================================
def mosaic_all_tiles_into_one_raster(src_names, dest_name):
    '''
    Description:
        Executes gdalwarp on the supplied source names to generate a mosaiced
        destination named file.
    '''

    logger = logging.getLogger(__name__)

    cmd = ['gdalwarp', '-wm', '2048', '-multi',
           '-srcnodata', str(NO_DATA_VALUE),
           '-dstnodata', str(NO_DATA_VALUE)]
    cmd.extend(src_names)
    cmd.append(dest_name)

    cmd = ' '.join(cmd)

    output = ''
    try:
        logger.info("Executing [{0}]".format(cmd))
        output = execute_cmd(cmd)
    except Exception:
        logger.error("Failed to mosaic tiles")
        raise
    finally:
        if len(output) > 0:
            logger.info(output)


# ============================================================================
def warp_raster_to_match_landsat_data(src_name, dest_name,
                                      src_proj4, dest_proj4,
                                      x_pixel_size, y_pixel_size,
                                      min_x_utm, min_y_utm,
                                      max_x_utm, max_y_utm):
    '''
    Description:
        Executes gdalwarp on the supplied source name to generate a warped
        destination named file that matches the projection and image extents
        and pixel size specified, which were derived from Landsat data.
    '''

    logger = logging.getLogger(__name__)

    cmd = ['gdalwarp', '-wm', '2048', '-multi',
           '-tr', str(x_pixel_size), str(y_pixel_size),
           '-s_srs', "'" + src_proj4 + "'",
           '-t_srs', "'" + dest_proj4 + "'",
           '-of', 'GTiff',
           '-overwrite', '-te',
           str(min_x_utm), str(min_y_utm), str(max_x_utm), str(max_y_utm),
           '-srcnodata', str(NO_DATA_VALUE),
           '-dstnodata', str(NO_DATA_VALUE)]
    cmd.append(src_name)
    cmd.append(dest_name)

    cmd = ' '.join(cmd)

    output = ''
    try:
        logger.info("Executing [{0}]".format(cmd))
        output = execute_cmd(cmd)
    except Exception:
        logger.error("Failed during warping to match Landsat")
        raise
    finally:
        if len(output) > 0:
            logger.info(output)


'''
TODO TODO TODO - NEED TO PROCESS A COASTAL SCENE BECAUSE MAY NOT HAVE ASTER TILE DATA FOR SOME OF THE SCENE AND WILL NEED TO DO SOMETHING ELSE
TODO TODO TODO - NEED TO PROCESS A COASTAL SCENE BECAUSE MAY NOT HAVE ASTER TILE DATA FOR SOME OF THE SCENE AND WILL NEED TO DO SOMETHING ELSE
TODO TODO TODO - NEED TO PROCESS A COASTAL SCENE BECAUSE MAY NOT HAVE ASTER TILE DATA FOR SOME OF THE SCENE AND WILL NEED TO DO SOMETHING ELSE
TODO TODO TODO - NEED TO PROCESS A COASTAL SCENE BECAUSE MAY NOT HAVE ASTER TILE DATA FOR SOME OF THE SCENE AND WILL NEED TO DO SOMETHING ELSE
TODO TODO TODO - NEED TO PROCESS A COASTAL SCENE BECAUSE MAY NOT HAVE ASTER TILE DATA FOR SOME OF THE SCENE AND WILL NEED TO DO SOMETHING ELSE
TODO TODO TODO - NEED TO PROCESS A COASTAL SCENE BECAUSE MAY NOT HAVE ASTER TILE DATA FOR SOME OF THE SCENE AND WILL NEED TO DO SOMETHING ELSE
TODO TODO TODO - NEED TO PROCESS A COASTAL SCENE BECAUSE MAY NOT HAVE ASTER TILE DATA FOR SOME OF THE SCENE AND WILL NEED TO DO SOMETHING ELSE
TODO TODO TODO - NEED TO PROCESS A COASTAL SCENE BECAUSE MAY NOT HAVE ASTER TILE DATA FOR SOME OF THE SCENE AND WILL NEED TO DO SOMETHING ELSE
TODO TODO TODO - NEED TO PROCESS A COASTAL SCENE BECAUSE MAY NOT HAVE ASTER TILE DATA FOR SOME OF THE SCENE AND WILL NEED TO DO SOMETHING ELSE
'''
# ============================================================================
def process(args):

    logger = logging.getLogger(__name__)

    # Read the XML metadata
    espa_xml = metadata_api.parse(args.xml_filename, silence=True)
    # Grab the global metadata object
    gm = espa_xml.get_global_metadata()
    # Grab the bands metadata object
    bands = espa_xml.get_bands()

    x_pixel_size = -1
    y_pixel_size = -1
    dest_proj4 = ''
    toa_bt_filename = ''
    # Find the TOA BT band to extract information from
    for band in bands.band:
        if band.product == "toa_bt":
            # Get the output pixel size
            x_pixel_size = band.pixel_size.x
            y_pixel_size = band.pixel_size.y

            toa_bt_filename = band.get_file_name()

            # Get the output proj4 string
            dest_proj4 = get_original_projection(toa_bt_filename)

            # Found the one we need to stop
            break

    # Error if we didn't find the TOA BT band in the data
    if x_pixel_size == -1:
        raise Exception("Failed to find the TOA BT band in the input data")

    # Determine the bounding geographic coordinates for the ASTER tiles we
    # will need
    east = math.ceil(gm.bounding_coordinates.east)
    west = math.floor(gm.bounding_coordinates.west)
    north = math.ceil(gm.bounding_coordinates.north)
    south = math.floor(gm.bounding_coordinates.south)

    # Determine the UTM projection corner points
    min_x_utm = 0
    max_x_utm = 0
    min_y_utm = 0
    max_y_utm = 0
    for cp in gm.projection_information.corner_point:
        if cp.location == 'UL':
            min_x_utm = cp.x
            max_y_utm = cp.y
        if cp.location == 'LR':
            max_x_utm = cp.x
            min_y_utm = cp.y

    # Adjust the UTM coordinates for image extents becuse they are in center
    # of pixel, and we need to supply the warpping with actual extents
    min_x_utm = min_x_utm - x_pixel_size * 0.5
    max_x_utm = max_x_utm + x_pixel_size * 0.5
    min_y_utm = min_y_utm - y_pixel_size * 0.5
    max_y_utm = max_y_utm + y_pixel_size * 0.5

    # Register all the gdal drivers and choose the GeoTiff for our temp output
    gdal.AllRegister()
    driver = gdal.GetDriverByName('GTiff')

    # The ASTER data is in geographic projection so specify that here
    ds_srs = osr.SpatialReference()
    ds_srs.ImportFromEPSG(4326)

    # Process through the lattitude and longitude ASTER tiles
    # - Download them
    # - Extract the Emissivity bands 13 and 14 as well as the NDVI
    # - Generate the B6 EMIS Landsat from the 13 and 14 band data
    emis_filenames = list()
    for lat in xrange (int(south), int(north)+1):
        for lon in xrange (int(west), int(east)+1):
            # Build the base filename using the positive or negative format
            filename = ''
            if lon < 0:
                filename = FILE_N_FORMAT.format(lat, lon)
            else:
                filename = FILE_P_FORMAT.format(lat, lon)

            # Build the HDF5 filename for the tile
            h5_file_path = ''.join([filename, '.h5'])

            # Build the URL and download the tile
            url_path = ''.join([SERVER, SERVER_PATH, h5_file_path])
            http_transfer_file(url_path, h5_file_path)

            # ----------------------------------------------------------------
            # Build the output tile names
            b6_emis_file_path = ''.join([filename, '_b6_emis.tif'])

            # Add the tile names to the list for mosaic building and warping
            emis_filenames.append(b6_emis_file_path)

            # Get the sub-dataset name for Emissivity-Mean
            emis_ds_name = ''.join(['HDF5:"', h5_file_path,
                                    '"://Emissivity/Mean'])
            # Get the sub-dataset name for Lattitude
            lat_ds_name = ''.join(['HDF5:"', h5_file_path,
                                   '"://Geolocation/Latitude'])
            # Get the sub-dataset name for Longitude
            lon_ds_name = ''.join(['HDF5:"', h5_file_path,
                                   '"://Geolocation/Longitude'])

            logger.debug(emis_ds_name)
            logger.debug(lat_ds_name)
            logger.debug(lon_ds_name)

            # ----------------------------------------------------------------
            try:
                (ds_13_band_data, ds_14_band_data, ds_lat_band_data,
                 ds_lon_band_data, x_dim, y_dim) = \
                    extract_aster_data(emis_ds_name, lat_ds_name, lon_ds_name)
            except Exception:
                logger.exception("Extracting ASTER data from tile")
                raise

            # Remove the HDF5 tile since we have extracted all the info we
            # need from it
            if os.path.exists(h5_file_path):
                os.unlink(h5_file_path)

            # ----------------------------------------------------------------
            # Determine the minimum and maximum latitude and longitude
            x_min = ds_lon_band_data.min()
            x_max = ds_lon_band_data.max()
            y_min = ds_lat_band_data.min()
            y_max = ds_lat_band_data.max()

            # Determine the resolution of the ASTER data
            x_res = (x_max - x_min) / float(x_dim)
            y_res = (y_max - y_min) / float(y_dim)

            # Build the geo transform to apply to the raster tile
            geo_transform = [x_min, x_res, 0, y_max, 0, -y_res]

            # ----------------------------------------------------------------
            # Create a mask of the band 13 data and then scale it
            # Original matlab code is <= 0
            ds_13_masked = np.ma.masked_where(ds_13_band_data <= 0,
                                              ds_13_band_data)
            ds_13_masked = ds_13_masked * 0.001
            del (ds_13_band_data)

            # Create a mask of the band 14 data and then scale it
            # Original matlab code is <= 0
            ds_14_masked = np.ma.masked_where(ds_14_band_data <= 0,
                                              ds_14_band_data)
            ds_14_masked = ds_14_masked * 0.001
            del (ds_14_band_data)

            # Create the Landsat B6 EMIS data
            ls_emis_b6 = 0.44 * ds_13_masked + 0.4 * ds_14_masked + 0.156

            # ----------------------------------------------------------------
            # Create the B6 EMIS raster output tile
            try:
                generate_raster_file(driver, b6_emis_file_path, x_dim, y_dim,
                                     ls_emis_b6.data, geo_transform,
                                     ds_srs.ExportToWkt(), NO_DATA_VALUE,
                                     gdal.GDT_Float32)
            except Exception:
                logger.exception("Generating B6 EMIS tile")
                raise

            # Memory cleanup
            del (ds_13_masked)
            del (ds_14_masked)

            del (ls_emis_b6)

    # Save the source proj4 string to use during warping
    src_proj4 = ds_srs.ExportToProj4()

    del (ds_srs)
    del (driver)

    # Define the temporary names
    temp_emis_name = 'temp_emis.tif'
    warped_emis_name = 'b6_emis.tif'

    # Mosaic the B6 EMIS tiles into the temp EMIS
    try:
        mosaic_all_tiles_into_one_raster(emis_filenames, temp_emis_name)
    except Exception:
        logger.exception("Mosaicing B6 EMIS tiles")
        raise

    # Cleanup the B6 EMIS tiles
    for emis_filename in emis_filenames:
        if os.path.exists(emis_filename):
            os.unlink(emis_filename)

    # Warp B6 EMIS to match the Landsat data
    try:
        warp_raster_to_match_landsat_data(temp_emis_name, warped_emis_name,
                                          src_proj4, dest_proj4,
                                          x_pixel_size, y_pixel_size,
                                          min_x_utm, min_y_utm,
                                          max_x_utm, max_y_utm)
    except Exception:
        logger.exception("Warping B6 EMIS to match Landsat data")
        raise

    # Cleanup the temp files
    if os.path.exists(temp_emis_name):
        os.unlink(temp_emis_name)

    # TODO TODO TODO
    # TODO TODO TODO
    # TODO TODO TODO
    # TODO TODO TODO
    # TODO TODO TODO
    # TODO TODO TODO
    # TODO TODO TODO
    # TODO TODO TODO
    # TODO TODO TODO
    # TODO TODO TODO
    # TODO TODO TODO
    # TODO TODO TODO
    # TODO TODO TODO
    # TODO TODO TODO
    # TODO TODO TODO


# ============================================================================
if __name__ == '__main__':
    '''
    Description:
      Generate B6 EMIS and ASTER NDVI from ASTER GED tiles for the specified
      Landsat scene.
    '''

    # Build the command line argument parser
    description = ("Retrieve ASTER data application")
    parser = ArgumentParser(description=description)

    # ---- Add parameters ----
    # Required parameters
    parser.add_argument('--xml',
                        action='store', dest='xml_filename', required=True,
                        help="The XML metadata file to use")

    # Parse the command line arguments
    args = parser.parse_args()

    # Configure logging
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:%(funcName)s'
                                ' -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        filename='aster.log',
                        level=logging.DEBUG)

    logger = logging.getLogger(__name__)

    try:
        # Call the main processing routine
        process(args)
    except Exception, e:
        if hasattr(e, 'output'):
            logger.error("Output [%s]" % e.output)
        logger.exception("Processing failed")
        sys.exit(EXIT_FAILURE)

    sys.exit(EXIT_SUCCESS)
