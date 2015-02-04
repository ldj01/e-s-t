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
def process_lst(args):
    '''
    Description:
        Provides the glue code for generating LST products.
    '''

    # Get the logger
    logger = logging.getLogger(__name__)

    # ------------------------------------------------------------------------
    # Retrieve and initial processing of the required AUX data
    cmd = ['extract_aux_data.py', '--xml', args.xml_filename]
    cmd = ' '.join(cmd)
    output = ''
    try:
        logger.info("Calling [{0}]".format(cmd))
        output = execute_cmd(cmd)
    except Exception, e:
        logger.error("Failed processing lst_download_extract_aux_data.py")
        raise e
    finally:
        if len(output) > 0:
            logger.info(output)

    # Extract the product id from the xml filename and build some other
    # filenames
    product_id = os.path.splitext(args.xml_filename)[0]
    mtl_filename = '{0}_MTL.txt'.format(product_id)
    dem_filename = '{0}_dem.img'.format(product_id)
    emi_filename = '{0}_emissivity.img'.format(product_id)

    # ------------------------------------------------------------------------
    # TODO TODO TODO - Get the DEM myself or does ESPA do that???
    # TODO TODO TODO - Get the DEM myself or does ESPA do that???
    # TODO TODO TODO - Get the DEM myself or does ESPA do that???
    cmd = ['do_create_dem.py',
           '--mtl', mtl_filename,
           '--dem', dem_filename]
    cmd = ' '.join(cmd)
    output = ''
    try:
        logger.info("Calling [{0}]".format(cmd))
        output = execute_cmd(cmd)
    except Exception, e:
        logger.error("Failed creating DEM")
        raise e
    finally:
        if len(output) > 0:
            logger.info(output)

    # ------------------------------------------------------------------------
    # TODO TODO TODO - The emnissivity data needs to be retrieved here and
    # TODO TODO TODO - processed.  Processed to required input???
    # TODO TODO TODO - The emnissivity data needs to be retrieved here and
    # TODO TODO TODO - processed.  Processed to required input???
    # TODO TODO TODO - The emnissivity data needs to be retrieved here and
    # TODO TODO TODO - processed.  Processed to required input???
    # TODO TODO TODO - The emnissivity data needs to be retrieved here and
    # TODO TODO TODO - processed.  Processed to required input???
    # TODO TODO TODO - The emnissivity data needs to be retrieved here and
    # TODO TODO TODO - processed.  Processed to required input???
    # TODO TODO TODO - The emnissivity data needs to be retrieved here and
    # TODO TODO TODO - processed.  Processed to required input???
    # TODO TODO TODO - The emnissivity data needs to be retrieved here and
    # TODO TODO TODO - processed.  Processed to required input???
    # TODO TODO TODO - The emnissivity data needs to be retrieved here and
    # TODO TODO TODO - processed.  Processed to required input???
    # ------------------------------------------------------------------------

    # ------------------------------------------------------------------------
    # Call the scene based lst
    cmd = ['lst',
           '--xml', args.xml_filename,
           '--dem', dem_filename,
           '--emi', emi_filename,
           '--verbose']
    cmd = ' '.join(cmd)
    output = ''
    try:
        logger.info("Calling [{0}]".format(cmd))
        output = execute_cmd(cmd)
    except Exception, e:
        logger.error("Failed processing scene_based_lst")
        raise e
    finally:
        if len(output) > 0:
            logger.info(output)

    # ------------------------------------------------------------------------
    # TODO TODO TODO - Don't know if anything here yet, except maybe cleanup
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
        Gathers input parameters and performs the LST processing.
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

    # Setup the default logger format and level.  Log to STDOUT.
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:'
                                '%(funcName)s -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO,
                        stream=sys.stdout)

    # Get the logger
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
