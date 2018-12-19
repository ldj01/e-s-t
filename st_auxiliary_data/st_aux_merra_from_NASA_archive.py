#! /usr/bin/env python3

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
import shutil
import sys
import logging
from datetime import datetime, timedelta

from espa.auxiliary.aux_utilities import System, SUCCESS, ERROR
from espa import AuxConfig
from espa import HttpSession
import cx_Oracle

# Global variables
MERRA_VARIABLES = ['H', 'T', 'QV']
logger = logging.getLogger(__name__)

def get_filename(cfg, year, month, day):
    """ Determines the base filename to process based on the
        date provided.

        Parameters:
            cfg: Configuration information
            year: Year for desired file
            month: Month for desired file
            day: Day for desired file

        Returns:
            Filename determined by the given date information
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

    name_format = cfg.get('nasa_filename_format')

    return name_format.format(stream, year, month, day)


def create_archive_file(source_filename, dest_filename):
    """ Extract the specified information into a new NetCDF4 file.  Only the
        variables needed by the Surface Temperature procedure is extracted,
        plus the indexes used to access those variables.

        Parameters:
            source_filename: Name of the source file to get the data from
            dest_filename: Name of the file to create

        Returns:
            SUCCESS/ERROR on success or failure of archive file creation
    """
    merra_variables_str = ",".join(MERRA_VARIABLES)
    cmd = ['nccopy -4 -V lon,lat,lev,time,' + merra_variables_str,
           source_filename, dest_filename]
    status = System.execute_cmd(' '.join(cmd))
    if status != SUCCESS:
        if os.path.exists(dest_filename):
            os.path.remove(dest_filename)
        return ERROR
    return SUCCESS


def extract_variables(cfg, source_filename, date):
    """ Extract the specified variables from the MERRA2 file and store
        that data in the archive

        Parameters:
            cfg: Configuration information
            source_filename: Name of the source file to get the data from
            date: Date for the file being processed

        Returns:
            SUCCESS/ERROR on success or failure of variable extraction
    """
    logger.debug("Processing [%s]", source_filename)

    # Figure out the filenames to create
    dload_dir = cfg.get('merra2_temp_directory')
    dest_filename = (cfg.get('merra2_filename_format')
                     .format(date.year, date.month, date.day, 'nc4'))
    dest_filename = os.path.join(dload_dir, dest_filename)
    logger.debug('Destination Filename [%s]', dest_filename)

    # Create archive file
    status = create_archive_file(source_filename, dest_filename)
    if status != SUCCESS:
        logger.warn('Processing for file %s to %s has failed. Continuing '
                    'processing.', source_filename, dest_filename)
        return ERROR

    # Archive the file
    arc_path = cfg.get('merra2_base_archive')
    arc_path = cfg.get('merra2_archive_format').format(arc_path, date.year)
    arc_file = os.path.join(arc_path, os.path.basename(dest_filename))
    logger.info('Archiving into [%s]', arc_path)
    shutil.move(dest_filename, arc_file)

    return SUCCESS


def archive_aux_data(args, cfg, archive, con, dbh):
    """ Iterates through dates and downloads, processes and archives the
        given date range

        Inputs:
            args: List of command line arguments
            cfg: Configuration information
            Archive: List of files that are in the current archive
            con: Connection to database
            dbh: Handle to interact with the database

        Returns:
            SUCCESS/ERROR on success or failure of creating new archive files
    """

    # Set up date variables
    s_date = args.start_date
    e_date = args.end_date

    # Establish a session and the login info
    auth = {}
    auth['username'] = os.getenv('NASA_EARTHDATA_USER')
    auth['password'] = os.getenv('NASA_EARTHDATA_PASS')
    if auth['username'] is None or auth['password'] is None:
        logger.error('Authentication for the NASA Earthdata site was not '
                     'provided.  Check the NASA_EARTHDATA_USER and '
                     'NASA_EARTHDATA_PASS environment variables.')
        return ERROR
    session = HttpSession(username=auth['username'], password=auth['password'])

    while s_date <= e_date:
        # Determine the directory to place the data and create it if it does
        # not exist
        arc_path = cfg.get('merra2_base_archive')
        arc_path = cfg.get('merra2_archive_format').format(arc_path, 
                s_date.year)
        System.create_directory(arc_path)

        # Determine archive filename for this date and if it needs processing
        arc_filename = cfg.get('merra2_filename_format').format(s_date.year,
                                                        s_date.month,s_date.day)
        arc_filename = os.path.join(arc_path, arc_filename)
        if not args.overwrite and arc_filename in archive:
            s_date += timedelta(days=1)
            continue

        # Determine NASA filename and url for this date
        nasa_filename = get_filename(cfg, s_date.year, s_date.month, s_date.day)
        if nasa_filename is None:
            logger.error('Could not generate the NASA filename.  Check the '
                         'config file for the nasa_filename_format variable')
            return ERROR

        data_url = cfg.get('nasa_url_format').format(s_date.year,
                                                     s_date.month,nasa_filename)
        logger.info('Retrieving %s', nasa_filename)

        # Add the temp download path to the filename
        dload_dir = cfg.get('merra2_temp_directory')
        System.create_directory(dload_dir)
        nasa_filename = os.path.join(dload_dir, nasa_filename)
        session.http_transfer_file(data_url, nasa_filename)

        logger.debug('Extracting variables')
        status = extract_variables(cfg, nasa_filename, s_date)
        if status != SUCCESS:
            # Processing failed, so continue to the next date
            s_date += timedelta(days=1)
            continue

        # Set up the variables that will be entered to the database
        values = {'coverage' : s_date}
        values['fname'] = os.path.basename(arc_filename)

        try:
            # Assume that if one file was in the archive, the group should also
            # be there
            existing = args.overwrite and arc_filename in archive
            if existing:
                # If we are overwriting and the file is already in the archive
                # then a db record already exists and must be updated
                logger.info('Updating existing database archive entry')
                dbh.execute('UPDATE MERRA2 SET FILE_NAME = :fname, '
                            'DATE_ENTERED = SYSDATE WHERE '
                            'EFFECTIVE_DATE = :coverage', values)
                if dbh.rowcount != 1:
                    logger.error('No rows updated')
                    raise cx_Oracle.DatabaseError
            else:
                logger.info('Creating new database archive entry')
                dbh.execute('INSERT INTO MERRA2 (EFFECTIVE_DATE, '
                            'FILE_NAME, DATE_ENTERED) VALUES (:coverage, '
                            ':fname, SYSDATE)', values)
            con.commit()
            logger.info('Processing successful for %s', s_date)
        except cx_Oracle.DatabaseError as e:
            logger.error('Database entry for %s has failed: %s'
                         'Exiting processing and cleaning the archive',
                         s_date, e)
            # If the file was not already in the archive, remove it
            if not existing:
                logger.error('Removing the archived file to match db records')
                os.remove(arc_filename)
            else:
                logger.error('This file previously existed in the local '
                             'archive. Double check that the archive and '
                             'database record for this file match.')

            return ERROR

        # Cleanup the download directory and increment date
        System.empty_directory(dload_dir)
        s_date += timedelta(days=1)
    # End s_date loop

    return SUCCESS


def main():
    """ Main routine which grabs the command-line arguments, determines
        which dates and times of data need to be processed, then processes the
        MERRA data that falls between the specified dates

        Returns:
            SUCCESS/ERROR on success or failure of MERRA-2 data retrieval and
            processing
    """

    # Get MERRA config information
    c = AuxConfig()
    cfg = c.get_config('merra2')
    start_date = datetime.strptime(cfg.get('merra2_start_date'), '%Y-%m-%d')

    # Gather command line arguments
    description = ('Downloads MERRA2 data and extracts the required parameters'
                   ' for Surface Temperature processing.  The parameters'
                   ' are then archived for later use.')
    args = System.get_command_line_arguments(description, start_date, False)

    # Setup logging
    System.setup_logging(args.debug)

    # Alert user if the overwrite flag is set
    if args.overwrite:
        logger.info('overwrite flag is set: Any existing archive data will '
                    'be overwritten by most current online archive data.')

    # Connect to the database
    try:
        con = cx_Oracle.connect(os.getenv('IAS_DB_COM'))
        dbh = con.cursor()
        logger.info('Connected to database successfully')
    except cx_Oracle.DatabaseError as e:
        logger.error('Could not connect to database: %s', e)
        sys.exit(ERROR)

    # Generate a list of the current archive files or create new archive
    archive = System.get_daily_archive_listing(args.start_date, args.end_date,
                                               cfg.get('merra2_base_archive'),
                                               cfg.get('merra2_archive_format'))

    status = archive_aux_data(args, cfg, archive, con, dbh)

    # Close the db connection and return
    dbh.close()
    con.close()
    sys.exit(status)

if __name__ == '__main__':
    main()
