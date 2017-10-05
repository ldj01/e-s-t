#! /usr/bin/env python

'''
    FILE: st_extract_auxiliary_narr_data.py

    PURPOSE: Performs setup and calls wgrib to extract NARR data from the
             binary auxiliary files.

    PROJECT: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    LICENSE: NASA Open Source Agreement 1.3
'''

import os
import sys
import logging
from argparse import ArgumentParser
from collections import namedtuple

from espa import Metadata
import st_utilities as util


PARMS_TO_EXTRACT = ['HGT', 'SPFH', 'TMP']
AUX_PATH_TEMPLATE = '{0:0>4}/{1:0>2}/{2:0>2}'
AUX_NAME_TEMPLATE = 'NARR_3D.{0}.{1:04}{2:02}{3:02}.{4:04}.{5}'
DATE_TEMPLATE = '{0:0>4}{1:0>2}{2:0>2}'


def retrieve_command_line_arguments():
    """Read arguments from the command line

    Returns:
        args <arguments>: The arguments read from the command line
    """

    parser = ArgumentParser(description='Retrieves and generates auxillary'
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
                            ('parameter', 'hdr', 'grb', 'output_dir'))


def build_aux_filenames(aux_path, parm, date, date_type):
    """Builds a filename set based on date and time

    Args:
        aux_path <str>: Path to base auxiliary (NARR) data
        parm <str>: NARR parameter to extract
        date <datetime>: Date and time of parameter to extract 
        date_type <str>: Type of date - time 0/1 (before/after scene center)

    Returns:
        <AuxFilenameSet>: The set of filenames
    """

    filename = AUX_NAME_TEMPLATE.format(parm,
                                        date.year,
                                        date.month,
                                        date.day,
                                        date.hour * 100,
                                        'hdr')

    path = AUX_PATH_TEMPLATE.format(date.year, date.month, date.day)

    hdr_path = os.path.join(aux_path, path, filename)

    return AuxFilenameSet(parameter=parm,
                          hdr=hdr_path,
                          grb=hdr_path.replace('.hdr', '.grb'),
                          output_dir='{0}_{1}'.format(parm, date_type))


def aux_filenames(aux_path, parms, t0_date, t1_date):
    """Builds t0 and t1 filename sets for each parameter

    Args:
        aux_path <str>: Path to base auxiliary (NARR) data
        parm <list[str]>: List of NARR parameters to extract
        t0_date <datetime>: Time 0 before scene center
        t1_date <datetime>: Time 1 after scene center

    Returns:
        list(<AuxFilenameSet>): A list of the filename sets
    """

    for parm in parms:
        yield build_aux_filenames(aux_path=aux_path, parm=parm,
                                  date=t0_date, date_type='t0')
        yield build_aux_filenames(aux_path=aux_path, parm=parm,
                                  date=t1_date, date_type='t1')


def extract_from_grib(aux_set):
    """Extracts the information from the grib file into a directory

    Args:
        aux_set <AuxFilenameSet>: Information needed to extract auxiliary data 
                                  from 1 file
    """

    logger = logging.getLogger(__name__)

    util.System.create_directory(aux_set.output_dir)

    with open(aux_set.hdr, 'r') as hdr_fd:
        for line in hdr_fd.readlines():
            parts = line.strip().split(':')
            record = parts[0]
            pressure = parts[6].split('=')[1]
            logger.debug('Processing Record {0}, Pressure {1}'
                         .format(record, pressure))

            path = os.path.join(aux_set.output_dir,
                                '.'.join([pressure, 'txt']))
            cmd = ['wgrib', aux_set.grb,
                   '-d', record,
                   '-text', '-nh', '-o', path]
            cmd = ' '.join(cmd)
            logger.info('wgrib command = [{}]'.format(cmd))

            # Extract the pressure data and raise any errors
            output = ''
            try:
                output = util.System.execute_cmd(cmd)
            except Exception:
                logger.error('Failed to unpack NARR Grib data')
                raise
            finally:
                if len(output) > 0:
                    logger.info(output)


def extract_narr_aux_data(espa_metadata, aux_path):
    """Extracts the required NARR data from the auxiliary archive

    Args:
        espa_metadata <espa.Metadata>: The metadata structure for the scene
        aux_path <str>: Path to base auxiliary (NARR) data
    """

    logger = logging.getLogger(__name__)

    (dummy, t0_date, t1_date) = util.NARR.dates(espa_metadata)

    logger.info('Before Date = {}'.format(str(t0_date)))
    logger.info(' After Date = {}'.format(str(t1_date)))

    for aux_set in aux_filenames(aux_path, PARMS_TO_EXTRACT,
                                 t0_date, t1_date):

        logger.info('Using {0}'.format(aux_set.hdr))
        logger.info('Using {0}'.format(aux_set.grb))

        # Verify that the files we need exist
        if (not os.path.exists(aux_set.hdr) or
                not os.path.exists(aux_set.grb)):
            raise Exception('Required ST AUX files are missing')

        extract_from_grib(aux_set)


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

    logger.info('*** Begin Extract Auxiliary NARR Data ***')

    # XML Metadata
    espa_metadata = Metadata()
    espa_metadata.parse(xml_filename=args.xml_filename)

    try:
        logger.info('Extracting ST AUX data')
        extract_narr_aux_data(espa_metadata, args.aux_path)

    except Exception:
        logger.exception('Failed processing auxiliary NARR data')
        raise

    logger.info('*** Extract Auxiliary NARR Data - Complete ***')


if __name__ == '__main__':
    main()
