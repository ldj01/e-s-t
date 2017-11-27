#! /usr/bin/env python


import os
import sys
import logging
import commands
import glob
import errno
from argparse import ArgumentParser


logger = None
# Name of the static NARR coordinates file
NARR_COORDINATES_PATH = ('../not-validated-prototype_st'
                         '/static_data/narr_coordinates.txt')


class LoggingFilter(logging.Filter):
    """Forces 'LSRD' to be provided in the 'subsystem' tag of the log format
       string
    """

    def filter(self, record):
        """Provide the string for the 'subsystem' tag"""

        record.subsystem = 'LSRD'

        return True


class ExceptionFormatter(logging.Formatter):
    """Modifies how exceptions are formatted
    """

    def formatException(self, exc_info):
        """Specifies how to format the exception text"""

        result = super(ExceptionFormatter,
                       self).formatException(exc_info)

        return repr(result)

    def format(self, record):
        """Specifies how to format the message text if it is an exception"""

        s = super(ExceptionFormatter, self).format(record)
        if record.exc_text:
            s = s.replace('\n', ' ')
            s = s.replace('\\n', ' ')

        return s


def setup_logging(args):
    """Configures the logging/reporting

    Args:
        args <args>: Command line arguments
    """

    global logger

    # Setup the logging level
    logging_level = logging.INFO
    if args.debug:
        logging_level = logging.DEBUG

    handler = logging.StreamHandler(sys.stdout)
    formatter = ExceptionFormatter(fmt=('%(asctime)s'
                                        ' %(subsystem)s'
                                        ' %(levelname)-8s'
                                        ' %(message)s'),
                                   datefmt='%Y-%m-%dT%H:%M:%S')

    handler.setFormatter(formatter)
    handler.addFilter(LoggingFilter())

    logger = logging.getLogger()
    logger.setLevel(logging_level)
    logger.addHandler(handler)


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


def get_arguments():
    parser = ArgumentParser(description='Extract Parameters')

    parser.add_argument('--debug',
                        action='store_true',
                        dest='debug',
                        required=False,
                        help='Debug logging')

    return parser.parse_args()


def process_parameter(filename, output_path, parameter, hour):

    output_base = '{0}/{1}-{2:0>2}'.format(output_path, parameter, hour)
    output_hdr = '.'.join([output_base, 'hdr'])
    output_grb = '.'.join([output_base, 'grb'])
    output_dat = '.'.join([output_base, 'dat'])
    output_good = '-good.'.join([output_base, 'txt'])
    output_fill = '-fill.'.join([output_base, 'txt'])

    cmds = [['wgrib', filename, '|', 'grep', parameter,
             '>', output_hdr],
            ['cat', output_hdr, '|', 'wgrib', filename, '-i', '-grib',
             '-o', output_grb],
            ['wgrib', output_grb, '|', 'grep', parameter, '>', output_hdr],
            ['wgrib', output_grb, '-d', '13', '-text', '-nh', '-o',
             output_dat]]

    for cmd in cmds:
        scmd = ' '.join(cmd)
        execute_cmd(scmd)

    coord_data = None
    with open(NARR_COORDINATES_PATH, 'r') as coords_fd:
        coord_data = [x.strip().split() for x in coords_fd]

    narr_data = None
    with open(output_dat, 'r') as narr_fd:
        narr_data = [x.strip() for x in narr_fd]

    with open(output_good, 'w') as good_fd:
        with open(output_fill, 'w') as fill_fd:
            for (l, r) in zip(coord_data, narr_data):
                l.append(r)
                value = float(r)
                if value != 9.999e+20:
                    good_fd.write(','.join(l))
                    good_fd.write('\n')
                else:
                    fill_fd.write(','.join(l))
                    fill_fd.write('\n')


def process_file(filename):
    print filename

    # Get the date information from the grib file
    parts = filename.split('.')
    year = int(parts[1][:4])
    month = int(parts[1][4:6])
    day = int(parts[1][6:8])
    hour = int(parts[1][8:])

    output_path = '{}.results'.format(filename)
    create_directory(output_path)

    for parameter in ['HGT', 'SPFH', 'TMP']:
        process_parameter(filename, output_path, parameter, hour)


def main():

    args = get_arguments()
    setup_logging(args)

    for merged_file in glob.glob('merged_AWIP32.*.3D'):
        process_file(merged_file)



if __name__ == '__main__':
    main()
