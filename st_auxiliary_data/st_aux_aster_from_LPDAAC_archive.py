#! /usr/bin/env python3

''' Retrieves archived ASTER files from the LP DAAC archive. It then extracts
    the variables required for ST processing (Emissivity/Mean, Emissivity/SDev,
    Geolocation/Latitude, Geolocation/Longitude, NDVI/Mean). After extraction,
    the file is compressed with h5repack, and archived to a location defined in
    level2_auxiliaries.conf.
'''

import os
import sys
import logging
from argparse import ArgumentParser

from espa.auxiliary.aux_utilities import System, SUCCESS, ERROR
from espa import AuxConfig
from espa import HttpSession

# Global variables
logger = logging.getLogger(__name__)
ASTER_VARIABLES = ['/Emissivity/Mean', '/Emissivity/SDev',
                   '/Geolocation/Latitude', '/Geolocation/Longitude',
                   '/NDVI/Mean']

def extract_and_archive(input_fname, output_fname, temp_path, archive_path):
    """ Takes an ASTER file, pulls out the ASTER_VARIABLES to a temp file,
        compresses it, and then archives the resultant file.

        Parameters:
            input_fname: Filename of the ASTER file retrieved from DAAC archive
            output_fname: Filename to store to in local archive
            temp_path: Path to the temporary download/processing directory
            archive_path: Path to the local archive

        Returns:
            SUCCESS/ERROR on successful or failed processing
    """
    # Define the name for the temp file and ensure it doesn't already exist
    temp_file = os.path.join(temp_path, 'aster_tmp.h5')
    if os.path.exists(temp_file):
        try:
            os.remove(temp_file)
        except Exception as e:
            logger.error('Error removing temporary file %s', temp_file)
            logger.error(e)
            return ERROR

    # Build and execute h5copy commands to extract the required variables
    for var in ASTER_VARIABLES:
        cmd = ['h5copy', '-i', os.path.join(temp_path, input_fname),
               '-o', temp_file, '-p', '-s', var, '-d', var]
        status = System.execute_cmd(' '.join(cmd))
        if status != SUCCESS:
            logger.error('Error returned from command [%s]', ' '.join(cmd))
            return ERROR

    # Compress the temp file into the archive and clean up the temp directory
    cmd = ['h5repack', '-i', temp_file, '-o',
           os.path.join(archive_path, output_fname), '-f', 'GZIP=1']
    status = System.execute_cmd(' '.join(cmd))
    if status != SUCCESS:
        logger.error('Error returned from command [%s]', ' '.join(cmd))
        return ERROR

    try:
        os.remove(temp_file)
        os.remove(os.path.join(temp_path, input_fname))
    except Exception as e:
        logger.error('Error removing temporary files from %s', temp_path)
        logger.error(e)
        return ERROR

    return SUCCESS


def archive_aux_data(cfg, tile_file):
    """ Creates a list of all ASTER filenames that need to be processed,
        downloads the DAAC archive file for each filename, and then processes
        that file.

        Parameters:
            cfg: Config information for the ASTER archive
            tile_file: File containing all ASTER filenames

        Returns:
            SUCCESS/ERROR on successful or failed processing
    """
    # Take out the temp dir and archive dir and ensure they exist
    temp_dir = cfg.get('aster_temp_directory')
    arc_dir = cfg.get('aster_archive_directory')
    try:
        System.create_directory(temp_dir)
        System.create_directory(arc_dir)
    except Exception as e:
        logger.error("Error setting up directories [%s] [%s]",temp_dir,arc_dir)
        logger.error(e)
        return ERROR

    # Get a list of files currently in the ASTER archive
    archive = System.directory_file_list(arc_dir)

    # Get a list of all required filenames
    with open(tile_file) as file_list:
        file_names = [line.strip() for line in file_list.readlines()]

    session = HttpSession(username=os.getenv('NASA_EARTHDATA_LOGIN'),
        password=os.getenv('NASA_EARTHDATA_PASSWORD'))

    # Go through each file name
    for fname in file_names:
        # Determine the new filename and ensure the archive doesn't contain it
        new_fname = fname.replace('h5', 'subset.h5')
        if os.path.join(arc_dir, new_fname) in archive:
            logger.error('The file %s is already in the ASTER archive. This '
                         'script needs to have an empty initial archive '
                         'directory to run.', fname)
            return ERROR

        logger.info('Processing file %s', fname)

        # Download the file from the DAAC archive
        success = session.http_transfer_file(cfg.get('aster_daac_url') + fname,
                                            os.path.join(temp_dir, fname))
        if not success:
            logger.error('Failed to download file [%s]', fname)
            logger.error('You may need to log in.  Check '
                '$NASA_EARTHDATA_LOGIN and '
                '$NASA_EARTHDATA_PASSWORD.')
            return ERROR

        # Extract variables and archive the new file
        status = extract_and_archive(fname, new_fname, temp_dir,
                                     cfg.get('aster_archive_directory'))
        if status != SUCCESS:
            logger.error('Failed to archive file [%s]', fname)
            return ERROR

    return SUCCESS

def main():
    """ Main routine which grabs the config information, grabs the ASTER files
        from the LP DAAC archive, extracts the variables required for ST
        processing, compresses the final file, and archives it to the directory
        specified in the config file.

        Returns:
            SUCCESS/ERROR on success or failure of ASTER data retrieval,
            processing, and archiving
    """

    # Grab the aster_ged_tile_list.txt file path from the command line
    parser = ArgumentParser('Downloads ASTER data from the LP DAAC archive, '
                            'then pulls out the data required for Level 2 '
                            'processing, and then compresses the file before '
                            'storing it in a local archive.')
    parser.add_argument('--aster-tile-list',
                        action='store', dest='aster_tile_list', required=True,
                        help=('File containing the list of all ASTER file '
                              'names.  The expected file is '
                              'aster_ged_tile_list.txt; however any file '
                              'containing the same information can be used.'))
    parser.add_argument('--debug',
                        action='store_true', dest='debug', default=False,
                        help='Provides debug logging.')

    args = parser.parse_args()

    # Setup logging
    System.setup_logging(args.debug)

    # Get ASTER config information
    c = AuxConfig()
    try:
        cfg = c.get_config('aster', True)
    except Exception as e:
        logger.error("Error getting aster configuration parameters")
        logger.error(e)

    sys.exit(archive_aux_data(cfg, args.aster_tile_list))

if __name__ == '__main__':
    main()
