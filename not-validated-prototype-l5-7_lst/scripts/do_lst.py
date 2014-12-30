#! /usr/bin/env python

'''
    FILE: do_lst.py

    PURPOSE: Calls the executables required to generate the LST products.

    PROJECT: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    LICENSE: NASA Open Source Agreement 1.3

    HISTORY:

    Date              Programmer               Reason
    ----------------  ------------------------ -------------------------------
    Jan/2015          Ron Dilley               Initial implementation
'''

import os
import sys
import errno
import re
import commands
import logging
import requests
from argparse import ArgumentParser

# espa-common objects and methods
from espa_constants import EXIT_FAILURE
from espa_constants import EXIT_SUCCESS
import metadata_api


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
def http_transfer_file(download_url, destination_file, headers=None):
    '''
    Description:
        Using http transfer a file from a source location to a destination
        file on the localhost.

        HTTP headers can be specified to modify how the download happens.
    '''

    logger = logging.getLogger(__name__)

    logger.info("Transfering {0}".format(download_url))

    req = requests.get(download_url, headers=headers)

    if not req.ok:
        logger.error("Transfer Failed - HTTP")
        req.raise_for_status()

    try:
        with open(destination_file, 'wb') as local_fd:
            local_fd.write(req.content)
    except:
        logger.error("Transfer Failed - HTTP")
        raise
    finally:
        req.close()
    logger.info("Transfer Complete - HTTP")


# ============================================================================
def process_lst(args):
    '''
    Description:
        Provides the glue code for generating LST products.
    '''

    # get the logger
    logger = logging.getLogger(__name__)

    # retrieve and initial processing of the required AUX data
    cmd = ['lst_download_extract_aux_data.py', '--xml', args.xml_filename]
    cmd = ' '.join(cmd)
    output = ''
    try:
        output = execute_cmd(cmd)
    except Exception, e:
        logger.error("Failed to unpack data")
        raise e
    finally:
        if len(output) > 0:
            logger.info(output)

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
        Performs gathers input parameters and performs the LST processing.
    '''

    # Create a command line arugment parser
    description = ("Retrieves and generates auxillary LST inputs, then"
                   " processes and calls other executables for LST generation")
    parser = ArgumentParser(description=description)

    # ---- Add parameters ----
    # Required parameters
    parser.add_argument('--xml',
                        action='store', dest='xml_filename', required=True,
                        help="The XML metadata file to use")

    # Parse the command line parameters
    args = parser.parse_args()

    # setup the default logger format and level. log to STDOUT.
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:'
                                '%(funcName)s -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO)

    # get the logger
    logger = logging.getLogger(__name__)

    if args.xml_filename == '':
        logger.fatal("No XML metadata filename provided.")
        logger.fatal("Error processing LST.  Processing will terminate.")
        sys.exit(EXIT_FAILURE)

    try:
        logger.info("Generating LST products")

        process_lst(args)

    except Exception, e:
        logger.exception("Error processing LST.  Processing will terminate.")
        sys.exit(EXIT_FAILURE)

    logger.info("Completion of LST processing")
    sys.exit(EXIT_SUCCESS)
