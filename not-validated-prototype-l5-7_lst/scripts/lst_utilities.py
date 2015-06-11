'''
    FILE: lst_utilities.py

    PURPOSE: Provide a library of routines to be used by LST python
             applications.

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
import commands
import requests
from time import sleep


# ============================================================================
class System:
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
class Web:
    '''
    Description:
        Provides methods for interfacing with web resources.
    '''

    # ------------------------------------------------------------------------
    @staticmethod
    def http_transfer_file(download_url, destination_file):
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
                req = session.get(url=download_url, timeout=300.0)

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

            except:
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
