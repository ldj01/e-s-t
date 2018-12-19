#! /usr/bin/env python3

'''
    PURPOSE: Retrieves archived NARR3D files from the CISL RDA for the dates
             requested.  Extracts the variables ST requires (HGT, TMP, SPFH)
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
import logging
import calendar
from datetime import datetime, timedelta

from espa.auxiliary.aux_utilities import System, SUCCESS, ERROR
from espa import AuxConfig
from espa import HttpSession
import cx_Oracle

# Global variables
NARR_VARIABLES = ['HGT', 'SPFH', 'TMP']
logger = logging.getLogger(__name__)

def gen_name_list(cfg, s_date, e_date, archive, overwrite):
    """ Determines all of the base filenames to download from the UCAR
        website.  This function will take the overwrite and archive parameters
        into account by not returning names of files that do not need to be
        reprocessed because they already exist in the archive and the overwrite
        flag is not set

        Parameters:
            cfg: Configuration information
            s_date: Starting date
            e_date: Ending date
            archive: List of files in current archive
            overwrite: Flag to determine if existing archive files should
                be overwritten or not

        Returns:
            Generator of base filenames to process

        Notes:
            Files typically contain 3 days.
            Special assumptions are coded for the end of the month; which
            may have 1, 2, 3, or 4 days depending on the month and year.
    """
    # Prepare archive and filename formats
    base_format = cfg.get('narr_base_archive')
    arc_format = cfg.get('narr_archive_format')
    ucar_file = cfg.get('ucar_filename_format')
    arc_file = cfg.get('narr_filename_format')

    c_date = s_date
    need_to_process = False

    while c_date <= e_date:
        fname = arc_format.format(base_format, c_date.year,
                                  c_date.month, c_date.day)
        fname = os.path.join(fname, arc_file.format('HGT', c_date.year,
                             c_date.month, c_date.day, c_date.hour*100, 'grb'))
        if overwrite or fname not in archive:
            need_to_process = True

        # Check that we are at the end of the hour range for the day, and that
        # the date is one that has a tar file - is a multiple of 3 except at the
        # end of the month
        if (need_to_process and c_date.hour == 21
           and ((c_date.day % 3 == 0 and c_date.day % 10 != 0)
           or c_date.day == calendar.monthrange(c_date.year, c_date.month)[1])):
            # Logic to determine if we are at the end of the month which
            # requires different logic for the start day
            if c_date.day==calendar.monthrange(c_date.year,c_date.month)[1]:
                start_day = 28
            else:
                start_day = c_date.day-2

            yield(ucar_file.format(c_date.year, c_date.month,
                                   start_day, c_date.day))
            need_to_process = False

        c_date += timedelta(hours=3)


def create_grib_hdr(grib_file, variable, hdr_name):
    """ Create an inventory/header file for the specified grib file

        Inputs:
            grib_file <str>: Name of the source grib file to get the header
            variable <str>: Variable to extract from the grib file
            hdr_name <str>: Name of the header file to create

        Returns:
            Status of the command executed to generate the file
    """
    cmd = ['wgrib', grib_file, '|', 'grep', variable, '>', hdr_name]
    status = System.execute_cmd(' '.join(cmd))
    return status


def create_grib_file(grib_file, hdr_name, grb_name):
    """ Create an inventory/header file for the specified grib file

        Inputs:
            grib_file: Name of the source grib file to get the data from
            hdr_name: Name of the header file for what to extract
            grb_name: Name of the grib file to create

        Returns:
            Status of the command executed to generate the file
    """
    cmd = ['cat', hdr_name, '|',
           'wgrib', grib_file, '-i', '-grib', '-o', grb_name]
    status = System.execute_cmd(' '.join(cmd))
    return status


def process_grib_for_variables(cfg, grib_file, archive, overwrite, con, dbh):
    """ Extract the specified variable from the grib file and archive it

        Inputs:
            cfg: Configuration information
            grib_file: Name of grib file to be processed into the archive
            archive: List of files in current archive
            overwrite: Flag to determine if existing archive files should
                be overwritten or not
            con: Connection to database
            dbh - Handle to interact with the database

        Returns:
            SUCCESS/ERROR on successful or failed variable processing
    """
    logger.debug("Processing [%s]", grib_file)

    # Get the date information from the grib file
    parts = grib_file.split('.') # Split filename from extension
    year = int(parts[1][:4])
    month = int(parts[1][4:6])
    day = int(parts[1][6:8])
    hour = int(parts[1][8:])
    date = datetime(year=year, month=month, day=day, hour=hour)

    # Create variables to keep track of new archive files and db command
    new_archive = []
    values = {'coverage' : date}

    # Generate temp and archive file paths
    dest_path = cfg.get('narr_base_archive')
    dest_path = cfg.get('narr_archive_format').format(dest_path,year,month,day)
    temp_dir = cfg.get('narr_temp_directory')

    for variable in NARR_VARIABLES:
        # Figure out the filenames to create
        hdr_name = (cfg.get('narr_filename_format')
                    .format(variable, year, month, day, hour*100, 'hdr'))
        grb_name = (cfg.get('narr_filename_format')
                    .format(variable, year, month, day, hour*100, 'grb'))
        # Check that this file is not already in the archive
        if not overwrite and os.path.join(dest_path, grb_name) in archive:
            logger.info('%s is already in the archive.  Skipping processing '
                        'for %s.',grb_name, date)
            return SUCCESS

        # Add the temp download dir path to the filename
        hdr_name = os.path.join(temp_dir, hdr_name)
        grb_name = os.path.join(temp_dir, grb_name)
        logger.debug('Hdr Name [%s]', hdr_name)
        logger.debug('Grb Name [%s]', grb_name)

        # Create inventory/header file to extract the variable data
        status = create_grib_hdr(grib_file, variable, hdr_name)
        if status != SUCCESS:
            logger.error('Error creating temporary inventory header file')
            return ERROR

        # Create grib file for the variable
        status = create_grib_file(grib_file, hdr_name, grb_name)
        if status != SUCCESS:
            logger.error('Error creating grib file')
            return ERROR

        # Create new inventory/header file for the variable
        status = create_grib_hdr(grb_name, variable, hdr_name)
        if status != SUCCESS:
            logger.error('Error creating header file')
            return ERROR

        # Keep track of new archive file info for later steps
        new_archive.append(os.path.basename(grb_name))
        new_archive.append(os.path.basename(hdr_name))
        values[variable + '_grib'] = os.path.basename(grb_name)
        # End variable loop

    # Create the archive directory if it doesn't exist
    System.create_directory(dest_path)

    # Archive the files
    for fname in new_archive:
        dest_file = fname
        logger.info('Archiving [%s] into [%s]', dest_file, dest_path)
        os.rename(os.path.join(temp_dir, dest_file),
                  os.path.join(dest_path, dest_file))

    # Update the database for the new files
    try:
        # Assume that if one file was in the archive, the group should also
        # be there
        existing = overwrite and os.path.join(dest_path,
                                              new_archive[0]) in archive
        if existing:
            # If we are overwriting and the file is already in the archive
            # then a db record already exists and must be updated
            logger.info('Updating existing database archive entry')
            dbh.execute('UPDATE NARR SET FILE_NAME_HEIGHT = :HGT_grib, '
                        'FILE_NAME_HUMIDITY = :SPFH_grib, FILE_NAME_TEMP = '
                        ':TMP_grib, DATE_ENTERED = SYSDATE WHERE '
                        'EFFECTIVE_DATE_TIME = :coverage', values)
            if dbh.rowcount != 1:
                logger.error('No rows updated')
                raise cx_Oracle.DatabaseError
        else:
            logger.info('Creating new database archive entry')
            dbh.execute('INSERT INTO NARR (EFFECTIVE_DATE_TIME, '
                        'FILE_NAME_HEIGHT, FILE_NAME_HUMIDITY, '
                        'FILE_NAME_TEMP, DATE_ENTERED) VALUES (:coverage, '
                        ':HGT_grib, :SPFH_grib, :TMP_grib, SYSDATE)', values)
        con.commit()
        logger.info('Processing successful for %s', date)
    except cx_Oracle.DatabaseError as e:
        logger.error('Database entry for %s has failed: %s'
                     'Exiting processing and cleaning the archive', date, e)
        # If the file was not already in the archive, remove it
        if not existing:
            logger.error('Removing the archived files to match db records')
            for fname in new_archive:
                os.remove(os.path.join(dest_path, fname))
        else:
            logger.error('These files previously existed in the local '
                         'archive. Double check that the archive and '
                         'database records for these files match.')

        return ERROR

    return SUCCESS

def archive_aux_data(s_date, e_date, cfg, overwrite, archive, con, dbh):
    """ Iterates through dates and downloads, processes, and archives the
        given dates

        Inputs:
            s_date: Start date for processing
            e_date: End date for processing
            cfg: Configuration information
            overwrite: Flag to determine if existing archive files should
            be overwritten or not
            archive: List of files in current archive
            con: Connection to database
            dbh: Handle to interact with the database

        Returns:
            SUCCESS/ERROR on successful or failed archive generation
    """
    # Create an HttpSession object for downloading then log in
    session = HttpSession()
    login_data = dict()
    login_data['action'] = 'login'
    login_data['email'] = os.getenv('UCAR_EMAIL')
    login_data['passwd'] = os.getenv('UCAR_PASSWORD')
    if login_data['email'] is None or login_data['passwd'] is None:
        logger.error('Login information was not found for the UCAR website.'
                     'Check the UCAR_EMAIL and UCAR_PASSWORD env variables.')
        return ERROR

    # Figure out the names of the files to retrieve
    names = gen_name_list(cfg, s_date, e_date, archive, overwrite)

    # Create the download directory if it doesn't exist
    temp_dir = cfg.get('narr_temp_directory')
    System.create_directory(temp_dir)

    # Establish a logged in session
    session.login(cfg.get('ucar_login_url'), login_data)

    for name in names:
        filename = '{0}.tar'.format(name)
        logger.info('Retrieving %s', filename)

        # Clear out the download directory for new batch of files
        System.empty_directory(temp_dir)

        year = name[7:11] # Parse the year out of the filename

        # Format the url and download the file
        data_url = cfg.get('ucar_url_format').format(year, filename)
        filename = os.path.join(temp_dir, filename)
        status = session.http_transfer_file(data_url, filename)

        if not status:
            logger.error('Unsuccessful download of %s to %s',data_url, filename)
            continue

        # Extract the tar'd data
        cmd = ['tar', '-xvf', filename, '-C', temp_dir]
        cmd = ' '.join(cmd)
        status = System.execute_cmd(cmd)
        if status != SUCCESS:
            logger.error('Failed to extract files from tarball')
            return ERROR

        # Get list of grib files
        grib_files = [x for x in System.directory_file_list(temp_dir)
                      if 'tar' not in x]

        # For each parameter we need
        for grib_file in grib_files:
            status = process_grib_for_variables(cfg, grib_file, archive,
                                                overwrite, con, dbh)
            if status != SUCCESS:
                logger.error('Unsuccessful processing of %s.', grib_file)
                return ERROR
    # End name loop

    return SUCCESS


def main():
    """ Main routine which grabs the command-line arguments, determines
        which dates and times of data need to be processed, then processes the
        NARR data that falls between the specified dates

        Returns:
            SUCCESS/ERROR on successful or failed processing

        Notes:
            1. Datetimes for which the archive already has NARR data will be
            skipped for download and processing unless the overwrite flag is
            set (default off).

            2. When checking if a file is in the existing archive, the HGT
            variable file will be checked as all 3 variables for a given date
            should be present when archived.
    """
    # Get NARR config information
    c = AuxConfig()
    cfg = c.get_config('narr')
    start_date = datetime.strptime(cfg.get('narr_start_date'), '%Y-%m-%d')

    # Gather command line arguments
    description = ('Downloads NARR data and extracts the required parameters'
                   ' for Surface Temperature processing.  The parameters'
                   ' are then archived for later use.  The NARR data is'
                   ' packaged into gzipped tar balls containing 1 to 4 days'
                   ' worth of data from the source site.  Because of that,'
                   ' all days contained in the package will always be'
                   ' processed.')
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

    # Gather a list of the current archive files
    archive = System.get_daily_archive_listing(args.start_date, args.end_date,
                                               cfg.get('narr_base_archive'),
                                               cfg.get('narr_archive_format'))

    # Download the files, archive them, and update the db archive records
    status = archive_aux_data(args.start_date, args.end_date, cfg,
                              args.overwrite, archive, con, dbh)

    # Close the db connection and return
    dbh.close()
    con.close()
    sys.exit(status)

if __name__ == '__main__':
    main()
