
'''
    PURPOSE: Provide a library of routines to be used by ST python
             applications.  Each routine is placed under a class in hopes of
             separating them into specific collections/groups.

    PROJECT: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    LICENSE: NASA Open Source Agreement 1.3
'''

import os
import logging
import errno
import commands
import datetime
from time import sleep
from cStringIO import StringIO
import requests
from osgeo import gdal, osr


class Version(object):
    '''
    Description:
        Provides methods for retrieving version information.
    '''

    version = '1.0.0'

    @staticmethod
    def version_number():
        '''
        Description:
            Returns the version number.
        '''

        return Version.version

    @staticmethod
    def version_text():
        '''
        Description:
            Returns the version information as a spelled out string.
        '''

        msg = 'Surface Temperature - Version {0}'.format(Version.version)
        return msg

    @staticmethod
    def app_version():
        '''
        Description:
            Returns the version information.
        '''

        version_text = 'st_{0}'.format(Version.version)
        return version_text


class System(object):
    '''
    Description:
        Provides methods for interfacing with the host server.
    '''

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

        logger.info('Executing [{0}]'.format(cmd))
        (status, output) = commands.getstatusoutput(cmd)

        if status < 0:
            message = 'Application terminated by signal [{0}]'.format(cmd)
            if len(output) > 0:
                message = ' Stdout/Stderr is: '.join([message, output])
            raise Exception(message)

        if status != 0:
            message = 'Application failed to execute [{0}]'.format(cmd)
            if len(output) > 0:
                message = ' Stdout/Stderr is: '.join([message, output])
            raise Exception(message)

        if os.WEXITSTATUS(status) != 0:
            message = ('Application [{0}] returned error code [{1}]'
                       .format(cmd, os.WEXITSTATUS(status)))
            if len(output) > 0:
                message = ' Stdout/Stderr is: '.join([message, output])
            raise Exception(message)

        return output

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

    @staticmethod
    def create_link(src_path, link_path):
        """Create the specified link with some error checking

        Args:
            src_path (str): The location where the link will point.
            link_path (str): The location where the link will reside.

        Raises:
            Exception()
        """

        # Create/Make sure the directory exists
        try:
            os.symlink(src_path, link_path)
        except OSError as ose:
            if (ose.errno == errno.EEXIST and os.path.islink(link_path) and
                    src_path == os.path.realpath(link_path)):
                pass
            else:
                raise


class NARR(object):
    """Provides common NARR data related methods
    """

    @staticmethod
    def dates(espa_metadata):
        """Determines the before(time_0), after(time_1), and aquisition dates

        Args:
            espa_metadata <espa.metadata>: The metadata for the data

        Returns:
            acquisition <datetime>: Scene center date and time
            time_0 <datetime>: NARR data datetime before scene center
            time_1 <datetime>: NARR data datetime after scene center
        """

        center_time = str(espa_metadata.xml_object
                          .global_metadata.scene_center_time)

        acq_date = str(espa_metadata.xml_object
                       .global_metadata.acquisition_date)

        # Join them while dropping the last two '<number>Z'
        date_time = '-'.join([acq_date, center_time[:-2]])

        acquisition = (datetime.datetime
                       .strptime(date_time, '%Y-%m-%d-%H:%M:%S.%f'))

        '''
        Determine the 3hr increments to use from the auxillary data
        We want the one before and after the scene acquisition time
        and convert back to formatted strings
        '''
        scene_hour = int(center_time[:2])
        time_0_hour = scene_hour - (scene_hour % 3)

        time_0 = datetime.datetime(acquisition.year,
                                   acquisition.month,
                                   acquisition.day,
                                   time_0_hour)

        time_1 = time_0 + datetime.timedelta(hours=3)

        # Round acquisition date to nearest minute.
        if acquisition.second >= 30:
            needed_seconds = datetime.timedelta(0, 60 - acquisition.second,
                -acquisition.microsecond)
            rounded_acquisition = acquisition + needed_seconds
        else:
            extra_seconds = datetime.timedelta(0, acquisition.second,
                acquisition.microsecond)
            rounded_acquisition = acquisition - extra_seconds

        return (rounded_acquisition, time_0, time_1)


class Web(object):
    '''
    Description:
        Provides methods for interfacing with web resources.
    '''

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
                    logger.error('HTTP - Transfer of [{0}] - FAILED'
                                 .format(download_url))
                    # The raise_for_status gets caught by this try's except
                    # block
                    req.raise_for_status()

                # Write the downloaded data to the destination file
                with open(destination_file, 'wb') as local_fd:
                    local_fd.write(req.content)

                # Break the looping
                done = True
                logger.info('HTTP - Transfer Complete')

            except Exception:
                logger.exception('HTTP - Transfer Issue')

                if req is not None:
                    status_code = req.status_code

                if status_code != requests.codes['not_found']:
                    if retry_attempt > 3:
                        logger.info('HTTP - Transfer Failed'
                                    ' - exceeded retry limit')
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


class Geo(object):
    '''
    Description:
        Provides methods for interfacing with geographic projections.
    '''

    @staticmethod
    def convert_imageXY_to_mapXY(image_x, image_y, transform):
        """Translate image coordinates into map coordinates"""

        map_x = (transform[0] +
                 image_x * transform[1] +
                 image_y * transform[2])
        map_y = (transform[3] +
                 image_x * transform[4] +
                 image_y * transform[5])

        return (map_x, map_y)

    @staticmethod
    def convert_mapXY_to_imageXY(map_x, map_y, transform):
        """Translate map coordinates into image coordinates"""

        # Convert the transform from image->map to map->image
        (success, inv_transform) = gdal.InvGeoTransform(transform)

        image_x = (inv_transform[0] +
                   map_x * inv_transform[1] +
                   map_y * inv_transform[2])
        image_y = (inv_transform[3] +
                   map_x * inv_transform[4] +
                   map_y * inv_transform[5])

        return (image_x, image_y)

    @staticmethod
    def get_proj4_projection_string(img_filename):
        '''
        Description:
            Determine the proj4 projection parameters for the specified image.

        Returns:
            proj4 - The proj4 projection string for the image.
        '''

        data_set = gdal.Open(img_filename)
        if data_set is None:
            raise RuntimeError('GDAL failed to open ({0})'
                               .format(img_filename))

        ds_srs = osr.SpatialReference()
        ds_srs.ImportFromWkt(data_set.GetProjection())

        proj4 = ds_srs.ExportToProj4()

        del ds_srs
        del data_set

        return proj4

    @staticmethod
    def update_envi_header(hdr_file_path, no_data_value):
        '''
        Description:
            Updates the specified ENVI header.  Especially the no data value,
            since it is not supported by the GDAL ENVI driver.
        '''

        hdr_text = StringIO()
        with open(hdr_file_path, 'r') as tmp_fd:
            while True:
                line = tmp_fd.readline()
                if not line:
                    break
                if (line.startswith('data ignore value') or
                        line.startswith('description')):
                    pass
                else:
                    hdr_text.write(line)

                if line.startswith('description'):
                    # This may be on multiple lines so read lines until
                    # we find the closing brace
                    if not line.strip().endswith('}'):
                        while 1:
                            next_line = tmp_fd.readline()
                            if (not next_line or
                                    next_line.strip().endswith('}')):
                                break
                    hdr_text.write('description ='
                                   ' {USGS-EROS-ESPA generated}\n')
                elif (line.startswith('data type') and
                      (no_data_value is not None)):
                    hdr_text.write('data ignore value = {0}\n'
                                   .format(no_data_value))

        # Do the actual replace here
        with open(hdr_file_path, 'w') as tmp_fd:
            tmp_fd.write(hdr_text.getvalue())

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
            del raster

        except Exception:
            raise

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
            logger.error('Failed to mosaic tiles')
            raise
        finally:
            if len(output) > 0:
                logger.info(output)
