'''
    FILE: lst_utilities.py

    PURPOSE: Provide a library of routines to be used by LST python
             applications.  Each routine is placed under a class in hopes of
             separating them into specific collections/groups.

    PROJECT: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    LICENSE: NASA Open Source Agreement 1.3

    HISTORY:

    Date              Programmer               Reason
    ----------------  ------------------------ -------------------------------
    June/2015         Ron Dilley               Initial implementation
'''

import os
import logging
import errno
import commands
import requests
from cStringIO import StringIO
from osgeo import gdal, osr
from time import sleep


# ============================================================================
class Version(object):
    '''
    Description:
        Provides methods for retrieving version information.
    '''

    version = '0.0.1'

    # ------------------------------------------------------------------------
    @staticmethod
    def version_text():
        msg = ('Landsat 5 and 7 - Land Surface Temperature - Version {0}'
               .format(Version.version))
        return msg

    # ------------------------------------------------------------------------
    @staticmethod
    def app_version():
        version_text = 'l5-7_lst_{0}'.format(Version.version)
        return version_text


# ============================================================================
class System(object):
    '''
    Description:
        Provides methods for interfacing with the host server.
    '''

    # ------------------------------------------------------------------------
    @staticmethod
    def execute_cmd(cmd):
        '''
        Description:
            Execute a command line and return the terminal output or raise an
            exception

        Returns:
            output - The stdout and/or stderr from the executed command.
        '''

        logger = logging.getLogger(__name__)

        output = ''

        logger.info("Executing [{0}]".format(cmd))
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

    # ------------------------------------------------------------------------
    @staticmethod
    def create_directory(directory):
        '''
        Description:
            Create the specified directory with some error checking.
        '''

        # Create/Make sure the directory exists
        try:
            os.makedirs(directory, mode=0755)
        except OSError as ose:
            if ose.errno == errno.EEXIST and os.path.isdir(directory):
                pass
            else:
                raise


# ============================================================================
class Web(object):
    '''
    Description:
        Provides methods for interfacing with web resources.
    '''

    # ------------------------------------------------------------------------
    @staticmethod
    def http_transfer_file(download_url, destination_file, headers=None):
        '''
        Description:
            Using http transfer a file from a source location to a destination
            file on the localhost.

        Returns:
            status_code - One of the following
                        - 200, requests.codes['ok']
                        - 404, requests.codes['not_found']:
                        - 503, requests.codes['service_unavailable']:

        Notes:
            If a 503 is returned, the logged exception should be reviewed to
            determine the real cause of the error.
        '''

        logger = logging.getLogger(__name__)

        logger.info(download_url)

        session = requests.Session()

        session.mount('http://', requests.adapters.HTTPAdapter(max_retries=3))
        session.mount('https://', requests.adapters.HTTPAdapter(max_retries=3))

        status_code = requests.codes['ok']
        retry_attempt = 0
        done = False
        while not done:
            status_code = requests.codes['ok']
            req = None
            try:
                req = session.get(url=download_url, timeout=300.0,
                                  headers=headers)

                if not req.ok:
                    logger.error("HTTP - Transfer of [{0}] - FAILED"
                                 .format(download_url))
                    # The raise_for_status gets caught by this try's except
                    # block
                    req.raise_for_status()

                # Write the downloaded data to the destination file
                with open(destination_file, 'wb') as local_fd:
                    local_fd.write(req.content)

                # Break the looping
                done = True
                logger.info("HTTP - Transfer Complete")

            except Exception:
                logger.exception("HTTP - Transfer Issue")

                if req is not None:
                    status_code = req.status_code

                if status_code != requests.codes['not_found']:
                    if retry_attempt > 3:
                        logger.info("HTTP - Transfer Failed"
                                    " - exceeded retry limit")
                        done = True
                    else:
                        retry_attempt += 1
                        sleep(int(1.5 * retry_attempt))
                else:
                    # Not Found - So break the looping because we are done
                    done = True

            finally:
                if req is not None:
                    req.close()

        return status_code


# ============================================================================
# TODO TODO TODO - This Warp class is not usable yet, it is only a very very
#                  early stage prototype.
class Warp(object):
    '''
    Description:
        Provides warping capabilities through GDAL's gdalwarp command line
        tool.
    '''

    supported_output_formats = ['ENVI', 'GTiff']
    x_pixel_min = 0.5
    x_pixel_max = 1000.0
    y_pixel_min = 0.5
    y_pixel_max = 1000.0

    # ------------------------------------------------------------------------
    def __init__(self):
        self.base_cmd = ['gdalwarp', '-wm', '2048', '-multi']

        self.x_pixel_size = None
        self.y_pixel_size = None

        self.source_proj4 = None
        self.target_proj4 = None

        self.overwrite_target = False

        self.target_image_extents = None

        self.output_format = None

        self.source_no_data_value = None
        self.target_no_data_value = None

        self.source_files = None
        self.target_file = None

    # ------------------------------------------------------------------------
    def set_output_pixel_size(x_pixel_size, y_pixel_size):
        if x_pixel_size < 0.5 or x_pixel_size > 1000.0:
            raise Exception('X pixel_size out of range [{0} - {1}]'
                            .format(self.x_pixel_min, self.x_pixel_max))

        if y_pixel_size < 0.5 or y_pixel_size > 1000.0:
            raise Exception('Y pixel_size out of range [{0} - {1}]'
                            .format(self.y_pixel_min, self.y_pixel_max))

        self.x_pixel_size = x_pixel_size
        self.y_pixel_size = y_pixel_size

    # ------------------------------------------------------------------------
    def set_source_proj4(source_proj4):
        self.source_proj4 = source_proj4

    # ------------------------------------------------------------------------
    def set_target_proj4(target_proj4):
        self.target_proj4 = target_proj4

    # ------------------------------------------------------------------------
    def set_source_no_data_value(no_data_value):
        self.source_no_data_value = float(no_data_value)

    # ------------------------------------------------------------------------
    def set_target_no_data_value(no_data_value):
        self.target_no_data_value = float(no_data_value)

    # ------------------------------------------------------------------------
    def overwrite_target(true_false):
        self.overwrite_target = true_false

    # ------------------------------------------------------------------------
    def set_output_format(format):

        if format not in self.supported_output_formats:
            raise NotImplementedError('Format [{0}] not supported'.
                                      format(format))

        self.output_format = format

    # ------------------------------------------------------------------------
    def set_source_files(source_files):
        self.source_files = source_files

    # ------------------------------------------------------------------------
    def set_target_file(target_file):
        self.target_file = target_file

    # ------------------------------------------------------------------------
    def execute(self):
        if self.source_files is None:
            raise Exception('Source file(s) not specified')

        if self.target_file is None:
            raise Exception('Target file not specified')

        # Add the base command
        cmd = [x for x in self.base_cmd]

        # Adde the pixel size command
        if self.x_pixel_size is not None and self.y_pixel_size is not None:
            cmd.extend(['-tr',
                        str(self.x_pixel_size),
                        str(self.y_pixel_size)])

        # TODO TODO TODO - Lots to get done here
        print ' '.join(cmd)


# ============================================================================
class Geo(object):
    '''
    Description:
        Provides methods for interfacing with web resources.
    '''

    # ------------------------------------------------------------------------
    @staticmethod
    def get_proj4_projection_string(img_filename):
        '''
        Description:
            Determine the proj4 projection parameters for the specified image.

        Returns:
            proj4 - The proj4 projection string for the image.
        '''

        ds = gdal.Open(img_filename)
        if ds is None:
            raise RuntimeError("GDAL failed to open (%s)" % img_filename)

        ds_srs = osr.SpatialReference()
        ds_srs.ImportFromWkt(ds.GetProjection())

        proj4 = ds_srs.ExportToProj4()

        del (ds_srs)
        del (ds)

        return proj4

    # ------------------------------------------------------------------------
    @staticmethod
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
                    sb.write('description = {USGS-EROS-ESPA generated}\n')
                elif (line.startswith('data type') and
                      (no_data_value is not None)):
                    sb.write('data ignore value = %s\n' % no_data_value)

        # Do the actual replace here
        with open(hdr_file_path, 'w') as tmp_fd:
            tmp_fd.write(sb.getvalue())

    # ------------------------------------------------------------------------
    @staticmethod
    def generate_raster_file(driver, filename, data, x_dim, y_dim,
                             geo_transform, proj_wkt,
                             no_data_value, data_type):
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
    @staticmethod
    def mosaic_tiles_into_one_raster(src_names, dest_name, no_data_value):
        '''
        Description:
            Executes gdalwarp on the supplied source names to generate a
            mosaic'ed destination named file.
        '''

        logger = logging.getLogger(__name__)

        cmd = ['gdalwarp', '-wm', '2048', '-multi',
               '-srcnodata', str(no_data_value),
               '-dstnodata', str(no_data_value)]
        cmd.extend(src_names)
        cmd.append(dest_name)

        cmd = ' '.join(cmd)

        output = ''
        try:
            output = System.execute_cmd(cmd)
        except Exception:
            logger.error("Failed to mosaic tiles")
            raise
        finally:
            if len(output) > 0:
                logger.info(output)
