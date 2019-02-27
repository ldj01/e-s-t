#! /usr/bin/env python3

'''
    FILE: st_extract_auxiliary_merra_data.py

    PURPOSE: Performs setup and extracts MERRA or GEOS5 data from the NetCDF
             auxiliary files.

    PROJECT: IAS/LPGS L1/L2 DPAS at the USGS EROS

    LICENSE: NASA Open Source Agreement 1.3
'''

import os
import sys
import logging
from argparse import ArgumentParser
from collections import namedtuple
import numpy as np
import netCDF4 as nc4

from espa import Metadata
import st_utilities as util

# H - Geopotential height
# QV - Specific humidity
# T - Air temperature
PARMS_TO_EXTRACT = ['H', 'QV', 'T']

PRESSURE_LAYERS = [1000, 975, 950, 925, 900, 875, 850, 825, 800,
                   775, 750, 725, 700, 650, 600, 550, 500, 450,
                   400, 350, 300, 250, 200, 150, 100, 70, 50, 40,
                   30, 20, 10, 7, 5, 4, 3, 2, 1, 0.7, 0.5, 0.4, 0.3, 0.1]

# Global template variables
AUX_PATH_TEMPLATE_GEOS5 = '{0:0>4}/{1:0>2}/{2:0>2}'
AUX_NAME_TEMPLATE_GEOS5 = 'fpit_{0:04}{1:02}{2:02}.{3:04}.{4}'

AUX_PATH_TEMPLATE_MERRA2 = '{0:0>4}/{1:0>2}'
AUX_NAME_TEMPLATE_MERRA2 = 'merra2_{0:04}{1:02}{2:02}.{3}'

def retrieve_command_line_arguments():
    """Read arguments from the command line

    Returns:
        args <arguments>: The arguments read from the command line
    """

    parser = ArgumentParser(description='Retrieves and generates auxiliary'
                                        ' ST inputs, then processes and'
                                        ' calls other executables for ST'
                                        ' generation')

    parser.add_argument('--version',
                        action='version',
                        version=util.Version.version_text())

    parser.add_argument('--xml',
                        action='store', dest='xml_filename',
                        required=False, default=None,
                        help='The XML metadata file to use')

    parser.add_argument('--aux_path',
                        action='store', dest='aux_path',
                        required=False, default=None,
                        help='Base auxiliary data directory to use')

    parser.add_argument('--reanalysis',
                        action='store', dest='reanalysis',
                        required=False, default='MERRA2',
                        choices=['NARR','MERRA2','GEOS5'],
                        help='Reanalysis source - NARR, MERRA2 or GEOS5')

    parser.add_argument('--debug',
                        action='store_true', dest='debug',
                        required=False, default=False,
                        help='Keep any debugging data')

    # Parse the command line parameters
    args = parser.parse_args()

    # Command line arguments are required so print the help if none were
    # provided
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)  # EXIT FAILURE

    if args.xml_filename is None:
        raise Exception('--xml must be specified on the command line')

    if args.xml_filename == '':
        raise Exception('No XML metadata filename provided.')

    if args.aux_path is None:
        raise Exception('--aux_path must be specified on the command line')

    if args.aux_path == '':
        raise Exception('No ST Auxiliary data directory provided.')

    return args


AuxFilenameSet = namedtuple('AuxFilenameSet',
                            ('parameter', 'nc4', 'hour', 'output_dir'))


def build_aux_filename(reanalysis, aux_path, parm, date, date_type):
    """Builds a filename set based on date and time

    Args:
        aux_path <str>: Path to base auxiliary (MERRA or GEOS5) data
        parm <str>: MERRA or GEOS5 parameter to extract
        date <datetime>: Date and time of parameter to extract 
        date_type <str>: Type of date - time 0/1 (before/after scene center)

    Returns:
        <AuxFilenameSet>: The set of filenames
    """

    if reanalysis == "GEOS5":
        filename = AUX_NAME_TEMPLATE_GEOS5.format(date.year,
                                                  date.month,
                                                  date.day,
                                                  date.hour * 100,
                                                  'nc4')
        path = AUX_PATH_TEMPLATE_GEOS5.format(date.year, date.month, date.day)

    elif reanalysis == "MERRA2":
        filename = AUX_NAME_TEMPLATE_MERRA2.format(date.year,
                                                   date.month,
                                                   date.day,
                                                   'nc4')
        path = AUX_PATH_TEMPLATE_MERRA2.format(date.year, date.month)

    nc4_path = os.path.join(aux_path, path, filename)
    return AuxFilenameSet(parameter=parm,
                          nc4=nc4_path,
                          hour=date.hour,
                          output_dir='{0}_{1}'.format(parm, date_type))

def aux_filenames(reanalysis, aux_path, parms, date, date_type):
    """Builds filename sets for each parameter

    Args:
        reanalysis <str>: Type of auxiliary data
        aux_path <str>: Path to base auxiliary (MERRA or GEOS5) data
        parms <list[str]>: List of MERRA or GEOS5 parameters to extract
        date <datetime>: Time to extract
        date_type <str>: Type of date - time 0/1 (before/after scene center)

    Returns:
        list(<AuxFilenameSet>): A list of the filename sets
    """

    for parm in parms:
        yield build_aux_filename(reanalysis, aux_path=aux_path, parm=parm,
                                  date=date, date_type=date_type)


def extract_from_netcdf(aux_set):
    """Extracts the information from the NetCDF file into a directory

    Args:
        aux_set <AuxFilenameSet>: Information needed to extract auxiliary data 
                                  from 1 file
    """

    logger = logging.getLogger(__name__)

    util.System.create_directory(aux_set.output_dir)

    # Read the NetCDF file.
    f = nc4.Dataset(aux_set.nc4, 'r')

    # Find the index for the hour. 
    # Datasets are indexed by [time][lev][lat][lon]
    hour_increment = int(24 / f.variables[aux_set.parameter].shape[0])
    hour_index = int(aux_set.hour / hour_increment)

    pressure_layers = f.variables['lev'][:]
    for pressure_index, pressure_layer in enumerate(PRESSURE_LAYERS):

        # Make sure the pressure layer that was read from the file matches
        # the value we're expecting
        if abs(pressure_layer - pressure_layers[pressure_index]) > 1.0e-06:
                raise Exception('Unexpected pressure layer {0} found in file'.
                        format(pressure_layers[pressure_index]))

        # Build the output file.
        output_filename = os.path.join(aux_set.output_dir, str(pressure_layer))

        # Write the data for the pressure layer, time, and parameter. 
        latlon = f.variables[aux_set.parameter][hour_index, pressure_index]

        if type(latlon) is np.ndarray: # no missing values, not masked array
            np.save(output_filename, latlon)
        else: # masked array
            np.save(output_filename, latlon.filled(9.999e+20))

def extract_merra_aux_data(espa_metadata, aux_path, reanalysis):
    """Extracts the required MERRA or GEOS5 data from the auxiliary archive

    Args:
        espa_metadata <espa.Metadata>: The metadata structure for the scene
        aux_path <str>: Path to base auxiliary (MERRA or GEOS5) data
    """

    logger = logging.getLogger(__name__)

    (dummy, t0_date, t1_date) = util.REANALYSIS.dates(espa_metadata)

    logger.info('Before Date = {}'.format(str(t0_date)))
    logger.info(' After Date = {}'.format(str(t1_date)))

    last_filename = ""
    for (date, date_type) in zip([t0_date, t1_date], ['t0', 't1']):
        for aux_set in aux_filenames(reanalysis, aux_path, PARMS_TO_EXTRACT,
                                     date, date_type):

            # Avoid logging the same filename repeatedly.
            if last_filename != aux_set.nc4:
                logger.info('Using {0}'.format(aux_set.nc4))
                last_filename = aux_set.nc4

            # Verify that the files we need exist
            if (not os.path.exists(aux_set.nc4)):
                raise Exception('Required ST AUX files are missing')

            extract_from_netcdf(aux_set)


def main():
    """Core processing for the application
    """

    # Command Line Arguments
    args = retrieve_command_line_arguments()

    # Setup the logging level
    log_level = logging.INFO
    if args.debug:
        log_level = logging.DEBUG

    # Setup the default logger format and level.  Log to STDOUT.
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:'
                                '%(funcName)s -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=log_level,
                        stream=sys.stdout)
    logger = logging.getLogger(__name__)

    logger.info('*** Begin Extract Auxiliary {0} Data ***'
               .format(args.reanalysis))

    # XML Metadata
    espa_metadata = Metadata()
    espa_metadata.parse(xml_filename=args.xml_filename)

    try:
        logger.info('Extracting ST AUX %s: data', args.reanalysis)
        extract_merra_aux_data(espa_metadata, args.aux_path, args.reanalysis)

    except Exception:
        logger.exception('Failed processing auxiliary {0} data'
                        .format(args.reanalysis))
        raise

    logger.info('*** Extract Auxiliary {0} Data - Complete ***'
               .format(args.reanalysis))


if __name__ == '__main__':
    main()
