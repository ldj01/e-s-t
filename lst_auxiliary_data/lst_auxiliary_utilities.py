#! /usr/bin/env python

'''
    PURPOSE: Provide a library of routines to be used by LST Auxiliary python
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
import requests
from time import sleep
from contextlib import closing
import json


# ============================================================================
class Version(object):
    '''
    Description:
        Provides methods for retrieving version information.
    '''

    version = '0.0.1'

    # ------------------------------------------------------------------------
    @staticmethod
    def version_number():
        '''
        Description:
            Returns the version number.
        '''

        return Version.version

    # ------------------------------------------------------------------------
    @staticmethod
    def version_text():
        '''
        Description:
            Returns the version information as a spelled out string.
        '''

        msg = ('Land Surface Temperature - Auxiliary Scripts - Version {0}'
               .format(Version.version))
        return msg


# ============================================================================
class Web(object):
    '''
    Description:
        Provides methods for interfacing with web resources.
    '''

    # ------------------------------------------------------------------------
    class Session(object):
        '''
        Description:
            Manages an http(s) session.
        '''

        # --------------------------------------------------------------------
        def __init__(self, max_retries=3, block_size=None, timeout=300.0):
            super(Web.Session, self).__init__()

            self.logger = logging.getLogger(__name__)

            self.session = requests.Session()

            self.timeout = timeout
            self.max_retries = max_retries
            self.status_code = requests.codes['ok']

            # Determine if we are streaming or not based on block_size
            self.block_size = block_size
            self.stream = False
            if self.block_size is not None:
                self.stream = True

            adapter = requests.adapters.HTTPAdapter(max_retries=max_retries)
            self.session.mount('http://', adapter)
            self.session.mount('https://', adapter)

        # --------------------------------------------------------------------
        def login(self, login_url, login_data):
            '''
            Description:
                Provides for establishing a logged in session.
            '''

            # Login to the site
            self.session.post(url=login_url, data=login_data)

        # --------------------------------------------------------------------
        def _get_file(self, download_url, destination_file, headers=None):
            '''
            Notes: Downloading this way will place the whole source file into
                   memory before dumping to the local file.
            '''

            with closing(self.session.get(url=download_url,
                                          timeout=self.timeout,
                                          headers=headers)) as req:

                self.status_code = req.status_code

                if not req.ok:
                    self.logger.error('HTTP - Transfer of [{0}] - FAILED'
                                      .format(download_url))
                    # The raise_for_status generates an exception to be caught
                    req.raise_for_status()

                # Write the downloaded data to the destination file
                with open(destination_file, 'wb') as local_fd:
                    local_fd.write(req.content)

        # --------------------------------------------------------------------
        def _stream_file(self, download_url, destination_file, headers=None):
            '''
            Notes: Downloading this way streams 'block_size' of data at a
                   time.
            '''

            retrieved_bytes = 0
            with closing(self.session.get(url=download_url,
                                          timeout=self.timeout,
                                          stream=True,
                                          headers=headers)) as req:
                self.status_code = req.status_code

                if not req.ok:
                    self.logger.error('HTTP - Transfer of [{0}] - FAILED'
                                      .format(download_url))
                    # The raise_for_status generates an exception to be caught
                    req.raise_for_status()

                file_size = int(req.headers['content-length'])

                # Set block size based on streaming
                if self.stream:
                    block_size = self.block_size
                else:
                    block_size = file_size

                # Write the downloaded data to the destination file
                with open(destination_file, 'wb') as local_fd:
                    for data_chunk in req.iter_content(block_size):
                        local_fd.write(data_chunk)
                        retrieved_bytes += len(data_chunk)

                if retrieved_bytes != file_size:
                    raise Exception('Transfer Failed - HTTP -'
                                    ' Retrieved {0} out of {1} bytes'
                                    .format(retrieved_bytes, file_size))

        # --------------------------------------------------------------------
        def http_transfer_file(self, download_url, destination_file):
            '''
            Description:
                Use http to transfer a file from a source location to a
                destination file on the localhost.
            Returns:
                status_code - One of the following
                            - 200, requests.codes['ok']
                            - 404, requests.codes['not_found']:
                            - 503, requests.codes['service_unavailable']:
            Notes:
                If a 503 is returned, the logged exception should be reviewed
                to determine the real cause of the error.
            '''

            self.logger.info(download_url)

            retry_attempt = 0
            done = False
            while not done:
                self.status_code = requests.codes['ok']
                try:

                    self._stream_file(download_url, destination_file)

                    self.logger.info("Transfer Complete - HTTP")
                    done = True

                except IOError:
                    self.logger.exception('HTTP - Transfer Issue')

                    if self.status_code not in (requests.codes['not_found'],
                                                requests.codes['forbidden']):

                        if retry_attempt > self.max_retries:
                            self.logger.info('HTTP - Transfer Failed'
                                             ' - exceeded retry limit')
                            done = True
                        else:
                            retry_attempt += 1
                            sleep(int(1.5 * retry_attempt))
                    else:
                        # Not Found - So break the looping because we are done
                        done = True

            return self.status_code

        # --------------------------------------------------------------------
        def get_lines_from_url(self, download_url):
            '''retrieve lines from a url'''

            data = []
            self.status_code = requests.codes['ok']

            with closing(self.session.get(url=download_url,
                                          timeout=self.timeout,
                                          stream=self.stream)) as req:
                self.status_code = req.status_code

                if not req.ok:
                    self.logger.error('HTTP - Transfer of [{0}] - FAILED'
                                      .format(download_url))
                    # The raise_for_status generates an exception to be caught
                    req.raise_for_status()

                for line in req.iter_lines():
                    data.append(line)

            return data


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


class Config(object):
    '''Provides access to configurable attributes of the script

    Provides transparent access to settings from configuration
        1.Settings are specified as a json object stored in a file.
            read_config will be used to insert these into Config.config dict.
        2.Settings can be defined in the dictionary, Config.config, as
            key/value pairs.

    Beware: read_config will overwrite the contents of the configuration file
    '''

    config_directory = None
    config_filename = 'lst_auxiliary.config'
    config_path = None
    config = None  # Holds result of reading json object from file.

    @classmethod
    def read_config(cls):
        '''Reads configurable options from a file'''

        # The config file is located in the same place as the executable
        if cls.config_directory is None:
            cls.config_directory = os.path.dirname(os.path.abspath(__file__))

        # Add the filename to the path
        if cls.config_path is None:
            cls.config_path = os.path.join(cls.config_directory,
                                           cls.config_filename)

        with open(cls.config_path, 'r') as config_fd:
            lines = list()
            for line in config_fd:
                # Skip rudimentary comments
                if line.strip().startswith('#'):
                    continue

                lines.append(line)

            cls.config = json.loads(' '.join(lines))

        if cls.config is None:
            raise Exception('Failed loading configuration')

        return cls.config

    @classmethod
    def get(cls, attribute_path):
        '''Get the value of a configurable setting

            First it will try to read json data from the config file.
            If file does not exist then onfig will be from default values.
            If Key/value pair doesn't exist in JSON-like object stored in the
                file then default values will be used.
        '''

        logger = logging.getLogger(__name__)

        if cls.config is None:
            cls.config = cls.read_config()

        logger.debug('Searching Config For - {0}'.format(attribute_path))

        config = cls.config
        for attribute in attribute_path.split('.'):
            logger.debug('Searching Config For Attribute - {0}'.format(attribute))
            if attribute in config:
                config = config[attribute]
            else:
                raise IOError('Configuration Item - {0} - Not Found'
                              .format(attribute_path))

        logger.debug('Found Config - {0}'.format(config))
        return config
