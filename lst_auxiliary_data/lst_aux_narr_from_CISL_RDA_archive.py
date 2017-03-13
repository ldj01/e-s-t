#! /usr/bin/env python

'''
    PURPOSE: Retieves archived NARR3D files from the CISL RDA for the dates
             requested.  Extracts the variables LST requires (HGT, TMP, SPFH)
             and repackages them into our internal location and filenames.

    PROJECT: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    LICENSE: NASA Open Source Agreement 1.3

    NOTES:

          NCEP     - National Centers for Environmental Prediction
                     http://www.ncep.noaa.gov

          NARR     - NCEP North American Regional Reanalysis

          CISL RDA - Computational & Information Systems Lab
                     Research Data Archive http://rda.ucar.edu

          NCAR     - National Center for Atmospheric Research
                     http://ncar.ucar.edu

          UCAR     - University Corporation for Atmospheric Research
                     http://www2.ucar.edu
'''


import os
import sys
import shutil
import logging
import calendar
from argparse import ArgumentParser
from datetime import datetime, timedelta

from lst_aux_exception import AuxiliaryError
from lst_aux_version import VERSION_TEXT
from lst_aux_logging import LoggingFilter, ExceptionFormatter
import lst_aux_config as config
from lst_aux_http_session import HttpSession
from lst_aux_parameters import NARR_VARIABLES

from lst_auxiliary_utilities import System


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


def determine_name_list(cfg, s_date, e_date):
    """Determines all of the base filenames to process based on the dates
       provided

    Args:
        cfg <ConfigInfo>: Configuration information
        s_date <datetime>: Starting date
        e_date <datetime>: Ending date

    Returns:
        <list>: Yields a list of base filenames to process

    Notes:
        Files typically contain 3 days.
        Special assumptions are coded for the end of the month; which
        may have 1, 2, 3, or 4 days ... depending on the month and year.
    """

    name_format = cfg.ucar.data_name_format

    days_3 = timedelta(days=3)
    c_date = s_date

    while c_date <= e_date:
        if c_date.day == 28:
            days = calendar.monthrange(c_date.year, c_date.month)[1]
            yield(name_format.format(c_date.year, c_date.month,
                                     c_date.day, days))
            delta = timedelta(days=(days - 28 + 1))
            c_date += delta
        else:
            yield(name_format.format(c_date.year, c_date.month,
                                     c_date.day, c_date.day + 2))
            c_date += days_3


def create_grib_hdr(grib_file, variable, hdr_name):
    """Create an inventory/header file for the specified grib file

    Args:
        grib_file <str>: Name of the source grib file to get the header
        variable <str>: Variable to extract fro mthe grib file
        hdr_name <str>: Name of the header file to create
    """

    cmd = ['wgrib', grib_file, '|', 'grep', variable, '>', hdr_name]
    output = System.execute_cmd(' '.join(cmd))
    if output is not None and len(output) > 0:
        logger.debug(output)


def create_grib_file(grib_file, hdr_name, grb_name):
    """Create extract the header specified information into a new grib file

    Args:
        grib_file <str>: Name of the source grib file to get the data from
        hdr_name <str>: Name of the header file for what to extract
        grb_name <str>: Name of the grib file to create
    """

    cmd = ['cat', hdr_name, '|',
           'wgrib', grib_file, '-i', '-grib', '-o', grb_name]
    output = System.execute_cmd(' '.join(cmd))
    if output is not None and len(output) > 0:
        logger.debug('Grib Contents ({})'
                     .format(', '.join(output.split())))


def process_grib_for_variable(cfg, variable, grib_file):
    """Extract the specified variable from the grib file and archive it
    """

    logger.debug("Processing [{0}]".format(grib_file))

    # Get the date information from the grib file
    parts = grib_file.split('.')
    year = int(parts[1][:4])
    month = int(parts[1][4:6])
    day = int(parts[1][6:8])
    hour = int(parts[1][8:])

    # Figure out the filenames to create
    hdr_name = (cfg.archive_name_format
                .format(variable, year, month, day, hour*100, 'hdr'))
    grb_name = (cfg.archive_name_format
                .format(variable, year, month, day, hour*100, 'grb'))
    logger.debug('Hdr Name [{}]'.format(hdr_name))
    logger.debug('Grb Name [{}]'.format(grb_name))

    # Create inventory/header file to extract the variable data
    create_grib_hdr(grib_file, variable, hdr_name)

    # Create grib file for the variable
    create_grib_file(grib_file, hdr_name, grb_name)

    # Create new inventory/header file for the variable
    create_grib_hdr(grb_name, variable, hdr_name)

    # Determine the directory to place the data and create it if it does
    # not exist
    dest_path = (cfg.archive_directory_format
                 .format(cfg.base_archive_directory, year, month, day))
    System.create_directory(dest_path)

    # Archive the files
    logger.info('Archiving into [{0}]'.format(dest_path))
    # Grib
    dest_file = os.path.join(dest_path, grb_name)
    shutil.copyfile(grb_name, dest_file)
    # Header
    dest_file = os.path.join(dest_path, hdr_name)
    shutil.copyfile(hdr_name, dest_file)

    # Cleanup the working directory
    if os.path.exists(grb_name):
        os.unlink(grb_name)
    if os.path.exists(hdr_name):
        os.unlink(hdr_name)


def archive_aux_data(args, cfg):
    """Provides the main archive processing
    """

    (s_date, e_date) = determine_date_range(args)

    # Figure out the names of the files to retrieve
    names = [x for x in determine_name_list(cfg, s_date, e_date)]
    logger.debug('Searching for ({})'.format(', '.join(names)))

    # Establish a logged in session
    if not args.dev_skip_download:
        session = HttpSession(block_size=cfg.transfer_block_size)

    # Log in
    login_data = dict()
    login_data['action'] = cfg.ucar.action
    login_data['email'] = cfg.ucar.email
    login_data['passwd'] = cfg.ucar.passwd
    if not args.dev_skip_download:
        session.login(cfg.ucar.login_url, login_data)

    for name in names:
        filename = '{0}.tar'.format(name)
        logger.info('Retrieving {0}'.format(filename))

        year = name[7:11]

        data_url = cfg.ucar.data_url_format.format(year, filename)

        if not args.dev_skip_download:
            session.http_transfer_file(data_url, filename)

            # Extract the tar'd data
        cmd = ['tar', '-xvf', filename]
        cmd = ' '.join(cmd)
        output = System.execute_cmd(cmd)
        if output is not None and len(output) > 0:
            logger.debug('Extracted Grib Files ({})'.format(', '.join(output.split())))
        grib_files = output.split()

        # For each parameter we need
        for variable in NARR_VARIABLES:
            logger.debug('Processing Variable [{0}]'.format(variable))
            for grib_file in grib_files:
                process_grib_for_variable(cfg, variable, grib_file)

        # Cleanup - Extracted grib files
        for grib_file in grib_files:
            if os.path.exists(grib_file):
                os.unlink(grib_file)

        # Cleanup - The Tar ball
        if not args.dev_keep_download and os.path.exists(filename):
            os.unlink(filename)


def arg_date_type(datestring):
    '''Validates the date string to be a specified format'''

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
    description = ('Downloads NARR data and extracts the required parameters'
                   ' for Land Surface Temperature processing.  The parameters'
                   ' are then archived for later use.  The NARR data is'
                   ' packaged into gzipped tar balls containing 1 to 4 days'
                   ' worth of data from the source site.  Because of that,'
                   ' all days contained in the package will always be'
                   ' processed.')
    parser = ArgumentParser(description=description)

    parser.add_argument('--version',
                        action='version',
                        version=VERSION_TEXT)

    parser.add_argument('--debug',
                        action='store_true',
                        dest='debug',
                        default=False,
                        help='display error information')

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

    # Start date must start on a day based on a 3day per file pattern
    day = (s_date.day - 1) / 3 * 3 + 1
    delta = timedelta(days=(day - s_date.day))
    s_date += delta

    return (s_date, e_date)


def main():
    """Provides the setup and executaion of the processor for the application
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
