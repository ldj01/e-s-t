#! /usr/bin/env python3

'''
    PURPOSE: Retrieves archived GEOS5 files from the NASA GES DISC for the
             dates requested. Extracts the variables ST requires (H, T, QV)
             and repackages them into our internal location and filenames.

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
GEOS_VARIABLES = ['H', 'T', 'QV']
logger = logging.getLogger(__name__)


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
    geos_variables_str = ",".join(GEOS_VARIABLES)
    cmd = ['nccopy -4 -V lon,lat,lev,time,' + geos_variables_str,
           source_filename, dest_filename]
    status = System.execute_cmd(' '.join(cmd))
    if status != SUCCESS:
        if os.path.exists(dest_filename):
            os.path.remove(dest_filename)
        return ERROR
    return SUCCESS


def extract_variables(cfg, source_filename, date, hour):
    """ Extract the specified variables from the GEOS5 file and store
        that data in the archive

        Parameters:
            cfg: Configuration information
            source_filename: Name of the source file to get the data from
            date: Date for the file being processed
            hour: Hour for the file being processed

        Returns:
            SUCCESS/ERROR on success or failure of variable extraction
    """
    logger.debug("Processing [%s]", source_filename)

    # Figure out the filenames to create
    dload_dir = cfg.get('geos5_temp_directory')
    dest_filename = (cfg.get('geos5_filename_format')
                     .format(date.year, date.month, date.day, hour))
    dest_filename = os.path.join(dload_dir, dest_filename)
    logger.debug('Destination Filename [%s]', dest_filename)

    # Create archive file
    status = create_archive_file(source_filename, dest_filename)
    if status != SUCCESS:
        logger.warn('Processing for file %s to %s has failed. Continuing '
                    'processing.', source_filename, dest_filename)
        return ERROR

    # Archive the file
    arc_path = cfg.get('geos5_base_archive')
    arc_path = cfg.get('geos5_archive_format').format(arc_path, date.year,
                                                      date.month, date.day)
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

    # Establish a session. Note the geos data comes from a special "backdoor"
    # where the username/password is not needed.  If they are defined, they
    # will be used, but are not required.
    auth = {}
    auth['username'] = os.getenv('NASA_EARTHDATA_USER')
    auth['password'] = os.getenv('NASA_EARTHDATA_PASS')
    session = HttpSession(username=auth['username'], password=auth['password'])

    while s_date <= e_date:
        # Convert year, month, day to day of year for nasa url
        doy = s_date.timetuple().tm_yday

        # Determine the directory to place the data and create it if it does
        # not exist
        arc_path = cfg.get('geos5_base_archive')
        arc_path = cfg.get('geos5_archive_format').format(arc_path,
                                                          s_date.year,
                                                          s_date.month,
                                                          s_date.day)
        System.create_directory(arc_path)

        hour = 0
        dload_dir = None
        for i in range(0, 8):
            # GEOS files are one file per 3 hours
            hour = (i * 3)
            if(s_date.year == e_date.year and
               s_date.month == e_date.month and
               s_date.day == e_date.day):
                # If retrieving today's files, make sure we are getting files
                # that should actually be available by comparing hours
                now = datetime.utcnow()
                # Guarantee 24 hour with %H
                now_hour = now.strftime('%H')
                if(hour > int(now_hour)):
                    logger.info(
                        'Cannot retrieve hour file ({0}) for today from {1}'.
                        format(hour*100, now))
                    break

            # Set date for entering into database
            db_date = datetime(year=s_date.year, month=s_date.month,
                               day=s_date.day, hour=hour)
            # Multiply by 100 to get the 4 digit non-zero hours for the
            # filenames
            hour *= 100

            # Determine archive filename for this date and if it needs
            # processing
            arc_filename = cfg.get('geos5_filename_format').format(
                                   s_date.year, s_date.month, s_date.day,
                                   hour)
            arc_filename = os.path.join(arc_path, arc_filename)
            if not args.overwrite and arc_filename in archive:
                # File in archive already and we are not overwriting, so
                # continue on to the next file
                logger.info(
                    'File in archive already - {0}.'.format(arc_filename))
                continue

            # Determine NASA filename and url for this date
            name_format = cfg.get('nasa_filename_format')
            nasa_filename = name_format.format(s_date.year, s_date.month,
                                               s_date.day, hour)
            if nasa_filename is None:
                logger.error('Could not generate the NASA filename. Check the '
                             'config file for the nasa_filename_format '
                             'variable')
                return ERROR

            data_url = cfg.get('nasa_url_format').format(s_date.year, doy,
                                                         nasa_filename)
            logger.info('Retrieving %s', nasa_filename)

            # Add the temp download path to the filename
            dload_dir = cfg.get('geos5_temp_directory')
            System.create_directory(dload_dir)
            nasa_filename = os.path.join(dload_dir, nasa_filename)
            session.http_transfer_file(data_url, nasa_filename)

            logger.debug('Extracting variables')
            status = extract_variables(cfg, nasa_filename, s_date, hour)
            if status != SUCCESS:
                # Processing failed for current file, so continue on to the
                # next file
                logger.error('Processing failed for {0}.'.
                             format(nasa_filename))
                continue

            # Set up the variables that will be entered to the database
            values = {'coverage': db_date}
            values['fname'] = os.path.basename(arc_filename)

            try:
                # If file is already in archive, it should also have an entry
                # in the database table
                existing = args.overwrite and arc_filename in archive
                if existing:
                    # If we are overwriting and the file is already in the
                    # archive, a db record already exists and must be updated
                    logger.info('Updating existing database archive entry')
                    dbh.execute('UPDATE GEOS5 SET FILE_NAME = :fname, '
                                'DATE_ENTERED = SYSDATE WHERE '
                                'EFFECTIVE_DATE_TIME = :coverage', values)
                    if dbh.rowcount != 1:
                        logger.error('No rows updated')
                        raise cx_Oracle.DatabaseError
                else:
                    logger.info('Creating new database archive entry')
                    dbh.execute('INSERT INTO GEOS5 (EFFECTIVE_DATE_TIME, '
                                'FILE_NAME, DATE_ENTERED) VALUES (:coverage, '
                                ':fname, SYSDATE)', values)
                con.commit()
                logger.info('Processing successful for %s %s', s_date, hour)
            except cx_Oracle.DatabaseError as e:
                logger.error('Database entry for %s-%s has failed: %s'
                             'Exiting processing and cleaning the archive',
                             s_date, hour, e)
                # If the file already exists in the archive, remove it since
                # the database entry failed so that the archive and the
                # database match
                if arc_filename in archive:
                    logger.error('Removing the archived file to match db '
                                 'records')
                    os.remove(arc_filename)
                else:
                    logger.error('This file previously existed in the local '
                                 'archive. Double check that the archive and '
                                 'database record for this file match.')

                return ERROR

        # Cleanup the download directory and increment date
        if(dload_dir is not None and os.path.exists(dload_dir)):
            System.empty_directory(dload_dir)
        s_date += timedelta(days=1)

    # End s_date loop

    return SUCCESS


def main():
    """ Main routine which grabs the command-line arguments, determines
        which dates and times of data need to be processed, then processes the
        GEOS data that falls between the specified dates

        Returns:
            SUCCESS/ERROR on success or failure of GEOS5 data retrieval and
            processing
    """

    # Get GEOS config information
    c = AuxConfig()
    cfg = c.get_config('geos5')
    start_date = datetime.strptime(cfg.get('geos5_start_date'), '%Y-%m-%d')

    # Gather command line arguments
    description = ('Downloads GEOS5 data and extracts the required parameters'
                   ' for Surface Temperature processing.  The parameters'
                   ' are then archived for later use.')
    args = System.get_command_line_arguments(description, start_date, False)

    # Setup logging
    System.setup_logging(args.debug)

    # Alert user if the overwrite flag is set
    if args.overwrite:
        logger.info('Overwrite flag is set: Any existing archive data will '
                    'be overwritten by most current online archive data.')

    # Connect to the database
    try:
        print ("IAS_DB_COM: {0}-{1}".format(os.getenv('IAS_DB_COM'),os.getenv('ORACLE_SID')))
        con = cx_Oracle.connect(os.getenv('IAS_DB_COM'))
        dbh = con.cursor()
        logger.info('Connected to database successfully')
    except cx_Oracle.DatabaseError as e:
        logger.error('Could not connect to database: %s', e)
        sys.exit(ERROR)

    # Generate a list of the current archive files or create new archive
    archive = System.get_daily_archive_listing(args.start_date, args.end_date,
                                               cfg.get('geos5_base_archive'),
                                               cfg.get('geos5_archive_format'))

    status = archive_aux_data(args, cfg, archive, con, dbh)

    # Close the db connection and return
    dbh.close()
    con.close()
    sys.exit(status)

if __name__ == '__main__':
    main()
