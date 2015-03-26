#! /usr/bin/env python

import os
import sys
import glob
import logging
import math
import requests
from cStringIO import StringIO
from argparse import ArgumentParser
from osgeo import gdal, osr
import numpy as np

# espa-common objects and methods
from espa_constants import EXIT_FAILURE
from espa_constants import EXIT_SUCCESS
import metadata_api


'''
 http://e4ftl01.cr.usgs.gov/ASTT/AG100.003/2000.01.01/AG100.v003.44.-077.0001.h5
'''
SERVER = 'http://e4ftl01.cr.usgs.gov'
SERVER_PATH = '/ASTT/AG100.003/2000.01.01/'
FILE_N_FORMAT = 'AG100.v003.{0:02}.{1:04}.0001'
FILE_P_FORMAT = 'AG100.v003.{0:02}.{1:03}.0001'

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
# END - http_transfer_file


def update_envi_header(hdr_file_path, no_data_value):
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
def process(args):

    logger = logging.getLogger(__name__)

    # Read the XML metadata
    espa_xml = metadata_api.parse(args.xml_filename, silence=True)
    # Grab the global metadata object
    gm = espa_xml.get_global_metadata()

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

    # Register all the gdal drivers and choose the GeoTiff for our temp output
    gdal.AllRegister()
    driver = gdal.GetDriverByName('GTiff')

    # The ASTER data is in geographic projection so specify that here
    ds_srs = osr.SpatialReference()
    ds_srs.ImportFromEPSG(4326)

    # Specify the no data value for the ASTER data
    no_data_value = -9999

    # Process through the lattitude and longitude ASTER tiles
    # - Download them
    # - Extract the Emissivity bands 13 and 14 as well as the NDVI
    # - Generate the B6 EMIS Landsat from the 13 and 14 band data
    emis_filenames = list()
    ndvi_filenames = list()
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
            # TEMP COMMENT FOR SPEED http_transfer_file(url_path, filename)

            # Build the output tile names
            b6_emis_file_path = ''.join([filename, '_b6_emis.tif'])
            ndvi_file_path = ''.join([filename, '_ndvi.tif'])

            # Add the tile names to the list for mosaic building and warping
            emis_filenames.append(b6_emis_file_path)
            ndvi_filenames.append(ndvi_file_path)

            # Get the sub-dataset name for Emissivity-Mean
            emis_ds_name = ''.join(['HDF5:"', h5_file_path,
                                    '"://Emissivity/Mean'])
            # Get the sub-dataset name for NDVI-Mean
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

            try:
                # Open the Emissivity sub-dataset
                emis_sds = gdal.Open(emis_ds_name)
                if emis_ds is None:
                    raise RuntimeError("GDAL failed to open"
                                       " {0}".format(emis_ds_name))

                # Open the NDVI sub-dataset
                ndvi_sds = gdal.Open(ndvi_ds_name)
                if ndvi_ds is None:
                    raise RuntimeError("GDAL failed to open"
                                       " {0}".format(ndvi_ds_name))

                # Open the Latitude sub-dataset
                lat_sds = gdal.Open(lat_ds_name)
                if lat_ds is None:
                    raise RuntimeError("GDAL failed to open"
                                       " {0}".format(lat_ds_name))

                # Open the Longitude sub-dataset
                lon_sds = gdal.Open(lon_ds_name)
                if lon_ds is None:
                    raise RuntimeError("GDAL failed to open"
                                       " {0}".format(lon_ds_name))

                # The dimensions are the same for all the bands so just use
                # the values from the Emissivity dataset
                x_dim = emis_sds.RasterXSize
                y_dim = emis_sds.RasterYSize

                # Retrieve the band 13 data from the HDF5 input
                ds_13_band_data = emis_sds.GetRasterBand(4).ReadAsArray(0, 0,
                                                                        x_dim,
                                                                        y_dim)

                # Retrieve the band 14 data from the HDF5 input
                ds_14_band_data = emis_sds.GetRasterBand(5).ReadAsArray(0, 0,
                                                                        x_dim,
                                                                        y_dim)
                del (emis_sds)

                # Retrieve the NDVI data from the HDF5 input
                ds_ndvi_band_data = ndvi_sds.GetRasterBand(1).ReadAsArray(0, 0,
                                                                          x_dim,
                                                                          y_dim)
                del (ndvi_sds)

                # Retrieve the Latitude data from the HDF5 input
                ds_lat_band_data = lat_sds.GetRasterBand(1).ReadAsArray(0, 0,
                                                                        x_dim,
                                                                        y_dim)
                del (lat_sds)

                # Retrieve the Longitude data from the HDF5 input
                ds_lon_band_data = lon_sds.GetRasterBand(1).ReadAsArray(0, 0,
                                                                        x_dim,
                                                                        y_dim)
                del (lon_sds)
            except Exception:
                raise

            # Determine the minimum and maximum latitude and longitude
            x_min = ds_lon_band_data.min()
            x_max = ds_lon_band_data.max()
            y_min = ds_lat_band_data.min()
            y_max = ds_lat_band_data.max()

            # Determine the resolution of the ASTER data
            x_res = (x_max - x_min) / float(x_dim)
            y_res = (y_max - y_min) / float(y_dim)

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

            # Create a mask of the NDVI band data and then scale it
            ds_ndvi_masked = np.ma.masked_where(ds_ndvi_band_data
                                                == no_data_value,
                                                ds_ndvi_band_data)
            ds_ndvi_masked = ds_ndvi_masked * 0.01
            del (ds_ndvi_band_data)

            # Create the B6 EMIS data
            ls_emis_b6 = 0.44 * ds_13_masked + 0.4 * ds_14_masked + 0.156

            # Build the transform to apply to the data
            geo_transform = [x_min, x_res, 0, y_max, 0, -y_res]

            try:
                # Create the B6 EMIS raster output tile
                emis_raster = driver.Create(b6_emis_file_path, x_dim, y_dim,
                                            1, gdal.GDT_Float32)

                emis_raster.SetGeoTransform(geo_transform)
                emis_raster.SetProjection(ds_srs.ExportToWkt())
                emis_raster.GetRasterBand(1).WriteArray(ls_emis_b6.data)
                emis_raster.GetRasterBand(1).SetNoDataValue(no_data_value)
                emis_raster.FlushCache()

                del (emis_raster)

                # Create the NDVI raster output tile
                ndvi_raster = driver.Create(ndvi_file_path, x_dim, y_dim,
                                            1, gdal.GDT_Float32)

                ndvi_raster.SetGeoTransform(geo_transform)
                ndvi_raster.SetProjection(ds_srs.ExportToWkt())
                ndvi_raster.GetRasterBand(1).WriteArray(ds_ndvi_masked.data)
                ndvi_raster.GetRasterBand(1).SetNoDataValue(no_data_value)
                ndvi_raster.FlushCache()

                del (ndvi_raster)
            except Exception:
                raise

            # Memory cleanup
            del (ds_13_masked)
            del (ds_14_masked)

            del (ds_ndvi_masked)

            del (ls_emis_b6)

    del (ds_srs)
    del (driver)

    print ' '.join(emis_filenames)
    print ' '.join(ndvi_filenames)

# Mosaic the tiles into one
#gdalwarp -wm 2048 -multi -t_srs '+proj=utm +zone=17 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs' AG100*ndvi.tif rdd.tif

# Warp to match the Landsat scene
#gdalwarp -wm 2048 -multi -s_srs '+proj=longlat +datum=WGS84 +no_defs ' -t_srs '+proj=utm +zone=17 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs' -tr 30 30 -dstnodata -9999 -te 573585.000 4678185.000 814215.000 4890615.000 rdd.tif rdd2.tif


# ============================================================================
if __name__ == '__main__':
    '''
    Description:
      Read parameters from the command line and build a JSON dictionary from
      them.  Pass the JSON dictionary to the process routine.
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
