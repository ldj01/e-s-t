#! /usr/bin/env python

'''
    PURPOSE: Retrieves archived MERRA-2 files from the NASA GES DISC for the 
             dates requested.  Extracts the variables ST requires (H, T, QV)
             and repackages them into our internal location and filenames.

    PROJECT: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    LICENSE: NASA Open Source Agreement 1.3

    NOTES:

          GES      - NASA Goddard Earth Sciences 
                     https://sciences.gsfc.nasa.gov/earth/

          DISC     - Data and Information Services Center 
                     https://disc.gsfc.nasa.gov/
'''


import os
import sys
import shutil
import logging
import calendar
from argparse import ArgumentParser
from datetime import datetime, timedelta

from st_aux_exception import AuxiliaryError
from st_aux_version import VERSION_TEXT
from st_aux_logging import LoggingFilter, ExceptionFormatter
import st_aux_merra_config as config
from st_aux_http_session import SessionWithHeaderRedirection
from st_aux_parameters import MERRA_VARIABLES
from st_aux_utilities import System


logger = None



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

    # Disable annoying INFO messages from the requests module
    logging.getLogger('requests').setLevel(logging.WARNING)
    # Disable annoying INFO messages from the urllib3 module
    logging.getLogger('urllib3').setLevel(logging.WARNING)


def get_filename(cfg, year, month, day):
    """Determines all of the base filenames to process based on the dates
       provided

    Args:
        cfg <ConfigInfo>: Configuration information
        year <int>: Starting year
        month <int>: Starting month
        day <int>: Starting day of month

    Returns:
        <list>: Yields a list of base filenames to process
    """

    # Assemble the production stream and version number portion of the
    # filename.  The time ranges for the streams is documented in
    # https://gmao.gsfc.nasa.gov/pubs/docs/Bosilovich785.pdf.
    # The stream is the first digit and version is currently always 00,
    # so we will just combine them.
    if year <= 1991:
        stream = 100
    elif year <= 2000:
        stream = 200
    elif year <= 2010:
        stream = 300
    else:
        stream = 400

    name_format = cfg.nasa.data_name_format

    return name_format.format(stream, year, month, day)


def determine_name_list(cfg, s_date, e_date):
    """Determines all of the base filenames to process based on the dates
       provided

    Args:
        cfg <ConfigInfo>: Configuration information
        s_date <datetime>: Starting date
        e_date <datetime>: Ending date

    Returns:
        <list>: Yields a list of base filenames to process
    """

    days_1 = timedelta(days=1)
    c_date = s_date

    while c_date <= e_date:
        yield(get_filename(cfg, c_date.year, c_date.month, c_date.day))
        c_date += days_1


def create_archive_file(source_filename, dest_filename):
    """Extract the specified information into a new NetCDF4 file.  Only the
       variables needed by the Surface Temperature procedure is extracted,
       plus the indexes used to access those variables.

    Args:
        source_filename <str>: Name of the source file to get the data from
        dest_filename <str>: Name of the file to create
    """

    merra_variables_str = ",".join(MERRA_VARIABLES)
    cmd = ['nccopy -4 -V lon,lat,lev,time,' + merra_variables_str, 
           source_filename, dest_filename]
    output = System.execute_cmd(' '.join(cmd))
    if output is not None and len(output) > 0:
        logger.debug('Contents ({})'
                     .format(', '.join(output.split())))


def extract_variables(cfg, source_filename):
    """Extract the specified variables from the MERRA2 file and archive it

    Args:
        cfg <ConfigInfo>: Configuration information
        source_filename <str>: Name of the source file to get the data from
    """

    logger.debug("Processing [{0}]".format(source_filename))

    # Get the date information from the source file
    parts = source_filename.split('.')
    year = int(parts[2][:4])
    month = int(parts[2][4:6])
    day = int(parts[2][6:8])

    # Figure out the filenames to create
    dest_filename = (cfg.archive_name_format
                .format(year, month, day, 'nc4'))
    logger.debug('Destination Filename [{}]'.format(dest_filename))

    # Create archive file
    create_archive_file(source_filename, dest_filename)

    # Determine the directory to place the data and create it if it does
    # not exist
    dest_path = (cfg.archive_directory_format
                 .format(cfg.base_archive_directory, year, month, day))
    System.create_directory(dest_path)

    # Archive the file
    logger.info('Archiving into [{0}]'.format(dest_path))
    dest_full_filename = os.path.join(dest_path, dest_filename)
    shutil.copyfile(dest_filename, dest_full_filename)

    # Cleanup the working directory
    if os.path.exists(dest_filename):
        os.unlink(dest_filename)


def archive_aux_data(args, cfg):
    """Provides the main archive processing

    Args:
        args <args>: Command line arguments
        cfg <ConfigInfo>: Configuration information
    """

    (s_date, e_date) = determine_date_range(args)

    # Figure out the names of the files to retrieve
    names = [x for x in determine_name_list(cfg, s_date, e_date)]
    logger.debug('Searching for ({})'.format(', '.join(names)))

    # Establish a logged in session
    if not args.dev_skip_download:
        session = SessionWithHeaderRedirection(cfg.nasa.username,
                                               cfg.nasa.password, 
                                               cfg.nasa.login_url)

    for name in names:
        logger.info('Retrieving {0}'.format(name))

        year = name[27:31]
        month = name[31:33]

        data_url = cfg.nasa.data_url_format.format(year, month, name)

        if not args.dev_skip_download:
            session.http_transfer_file(data_url, name)

        logger.debug('Extracting variables')
        extract_variables(cfg, name)

        # Cleanup the original downloaded file 
        if not args.dev_keep_download and os.path.exists(name):
            os.unlink(name)


def arg_date_type(datestring):
    """Validates the date string to be a specified format

    Args:
        datestring <str>: String entered by user in a date field 
    """

    try:
        return datetime.strptime(datestring, '%Y%m%d').date()
    except ValueError:
        print('Dates must be the in the format: "YYYYMMDD"')
        raise


def get_command_line_arguments():
    """Setup and return the command line options

    Returns:
        <args>: Command line arguments
    """

    # Create a command line arugment parser
    description = ('Downloads MERRA2 data and extracts the required parameters'
                   ' for Surface Temperature processing.  The parameters'
                   ' are then archived for later use.')
    parser = ArgumentParser(description=description)

    parser.add_argument('--version',
                        action='version',
                        version=VERSION_TEXT)

    parser.add_argument('--debug',
                        action='store_true',
                        dest='debug',
                        default=False,
                        help='provides debug logging')

    parser.add_argument('--start-date',
                        action='store', dest='start_date',
                        metavar='YYYYMMDD', type=arg_date_type,
                        required=False,
                        help=('The start date YYYYMMDD(inclusive)'
                              ' if specifying a range'
                              '  Defaults to --end-date if not specified'))

    parser.add_argument('--end-date',
                        action='store', dest='end_date',
                        metavar='YYYYMMDD', type=arg_date_type,
                        required=False,
                        help=('The end date YYYYMMDD(inclusive)'
                              ' if specifying a range'))

    parser.add_argument('--date',
                        action='store', dest='date',
                        metavar='YYYYMMDD', type=arg_date_type,
                        required=False,
                        help='The date YYYYMMDD for a specific date')

    parser.add_argument('--dev-skip-download',
                        action='store_true', dest='dev_skip_download',
                        required=False,
                        help='Skip downloading')

    parser.add_argument('--dev-keep-download',
                        action='store_true', dest='dev_keep_download',
                        required=False,
                        help='Keep downloaded files')

    # Parse the command line parameters
    args = parser.parse_args()

    return args


def determine_date_range(args):
    """Determines the date range to use during the retrieval

    Args:
        <args>: Command line arguments
    """

    s_date = None
    e_date = None

    # Check if dates were given
    if args.date is not None:
        s_date = args.date
        e_date = args.date

    elif args.end_date is not None:
        e_date = args.end_date

        if args.start_date is not None:
            s_date = args.start_date

        else:
            s_date = e_date
    else:
        raise AuxiliaryError('Must supply either --date or --end-date')

    if e_date < s_date:
        raise AuxiliaryError('--end-date must be equal-to or after --start-date')

    return (s_date, e_date)


def main():
    """Provides the setup and execution of the processor for the application
    """

    args = get_command_line_arguments()

    setup_logging(args)

    cfg = config.get_config()

    archive_aux_data(args, cfg)


if __name__ == '__main__':
    try:
        main()
    except Exception:
        if logger is not None:
            logger.exception('EXCEPTION ')
