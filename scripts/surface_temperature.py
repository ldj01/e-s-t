#! /usr/bin/env python

'''
    PURPOSE: Determine which executable to run and then pass all arguments
             through to the appropriate script.

    PROJECT: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    LICENSE: NASA Open Source Agreement 1.3

    NOTES:
        This script does not have its own help message and will just return
        the help from underlying executables where appropriate.

        If this script has a required argument then only the usage for that
        argument will be shown if that argument is not included.

        All output from the underlying script will be given to the logger as
        an info message.
'''

import os
import sys
import logging
import argparse
import commands
from ConfigParser import ConfigParser


class ExecuteError(Exception):
    """Raised when command in execute_cmd returns with error"""

    def __init__(self, message, *args):
        self.message = message
        Exception.__init__(self, message, *args)


def execute_cmd(cmd):
    """Execute a command line and return the terminal output

    Args:
        cmd <str>: The command string to execute

    Returns:
        <str>: The stdout and/or stderr from the executed command

    Raises:
        ExecuteError(<str>)
    """

    (status, output) = commands.getstatusoutput(cmd)

    if status < 0:
        message = ('Application terminated by signal [{0}]'
                   .format(cmd))
        if len(output) > 0:
            message = ' Stdout/Stderr is: '.join([message, output])
        raise ExecuteError(message)

    if status != 0:
        message = 'Application failed to execute [{0}]'.format(cmd)
        if len(output) > 0:
            message = ' Stdout/Stderr is: '.join([message, output])
        raise ExecuteError(message)

    if os.WEXITSTATUS(status) != 0:
        message = ('Application [{0}] returned error code [{1}]'
                   .format(cmd, os.WEXITSTATUS(status)))
        if len(output) > 0:
            message = ' Stdout/Stderr is: '.join([message, output])
        raise ExecuteError(message)

    return output


def parse_cmd_line():
    """Will only parse --xml XML_FILENAME from cmdline.

    Precondition:
        '--xml FILENAME' exists in command line arguments

    Postcondition:
        returns xml_filename

    Note: Help is not included because the program will return the help from
          the underlying program.
    """

    # Try to parse out the XML so the application can be determined
    parse_xml = argparse.ArgumentParser(add_help=False)
    parse_xml.add_argument('--xml', action='store',
                           dest='xml_filename', required=True,
                           help='Input XML metadata file',
                           metavar='FILE')

    (temp, extra_args) = parse_xml.parse_known_args()

    return temp.xml_filename


def get_satellite_sensor_code(xml_filename):
    """Returns the satellite-sensor code if known
    """

    old_prefixes = ['LT4', 'LT5', 'LE7', 'LT8', 'LC8', 'LO8']
    collection_prefixes = ['LT04', 'LT05', 'LE07', 'LT08', 'LC08', 'LO08']

    base_name = os.path.basename(xml_filename)

    satellite_sensor_code = base_name[0:3]
    if satellite_sensor_code in old_prefixes:
        return satellite_sensor_code

    satellite_sensor_code = base_name[0:4]
    if satellite_sensor_code in collection_prefixes:
        return satellite_sensor_code

    raise Exception('Satellite-Sensor code ({0}) not understood'
                    .format(satellite_sensor_code))


def get_science_application_name(satellite_sensor_code):
    """Returns name of executable that needs to be called
    """

    available = ['LT4', 'LT5', 'LE7', 'LC8',
                 'LT04', 'LT05', 'LE07', 'LC08']

    if satellite_sensor_code in available:
        return 'st_generate_products.py'
    else:
        raise Exception('Satellite-Sensor code ({0}) not understood'
                        .format(satellite_sensor_code))


def main():
    """Determines executable, and calls it with all input arguments
    """

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

    xml_filename = parse_cmd_line()
    satellite_sensor_code = get_satellite_sensor_code(xml_filename)

    # Get the science application
    cmd = [get_science_application_name(satellite_sensor_code)]
    # Pass all arguments through to the science application
    cmd.extend(sys.argv[1:])

    # Convert the list to a string
    cmd = ' '.join(cmd)
    try:
        logger.info(' '.join(['EXECUTING SCIENCE APPLICATION:', cmd]))
        output = execute_cmd(cmd)

        if len(output) > 0:
            logger.info('\n{0}'.format(output))
    except ExecuteError:
        logger.exception('Error running {0}.'
                         'Processing will terminate.'
                         .format(os.path.basename(__file__)))
        raise  # Re-raise so exception message will be shown.


if __name__ == '__main__':
    main()
