#! /usr/bin/env python

'''
    FILE: lst_extract_tape6_results.py

    PURPOSE: Parses the tape6 and extracts the WAVELENGTH and TOTAL_RADIANCE
             values.  Writes the values to the requested output file.

    PROJECT: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    LICENSE: NASA Open Source Agreement 1.3

    HISTORY:

    Date              Programmer               Reason
    ----------------  ------------------------ -------------------------------
    Jan/2015          Ron Dilley               Initial implementation
'''

import os
import sys
import logging
from argparse import ArgumentParser
from cStringIO import StringIO


# espa-common objects and methods
from espa_constants import EXIT_FAILURE
from espa_constants import EXIT_SUCCESS


# ============================================================================
def extract_tpst(tape6_fd):
    '''
    Description:
        Extract the target pixel surface temperature value from the tape6_fd
        file and return it as a float value.
    '''

    target_pixel_surface_temp = None

    # Search for the TARGET-PIXEL (H2) SURFACE TEMPERATURE
    for line in tape6_fd:
        # Remove all whitespace and newlines
        line = ' '.join(line.strip().split())

        if line.startswith('TARGET-PIXEL (H2) SURFACE TEMPERATURE [K]'):
            target_pixel_surface_temp = float(list(reversed(line.split()))[0])
            break

    return target_pixel_surface_temp


# ============================================================================
def create_output(args, target_pixel_surface_temp, record_count, records):
    '''
    Description:
        Creates the output file.
    '''

    radiance_dat_filename = '{0}/lst_modtran.dat'.format(args.output_path)
    radiance_info_filename = '{0}/lst_modtran.info'.format(args.output_path)

    with open(radiance_dat_filename, "w") as radiance_dat_fd:
        radiance_dat_fd.write(records.getvalue())

    with open(radiance_info_filename, "w") as radiance_info_fd:
        radiance_info_fd.write("TARGET_PIXEL_SURFACE_TEMPERATURE {0}\n"
                               .format(target_pixel_surface_temp))
        radiance_info_fd.write("RADIANCE_RECORD_COUNT {0}\n"
                               .format(record_count))


# ============================================================================
def process_tape6_results(args):
    '''
    Description:
        Parses the tape6 and extracts the WAVELENGTH and TOTAL_RADIANCE
        values.  Writes the values to the requested output file.
    '''

    logger = logging.getLogger(__name__)

    tape6_filename = '{0}/tape6'.format(args.input_path)

    if not os.path.isfile(tape6_filename):
        raise Exception("Missing tape6 file in input directory")

    target_pixel_surface_temp = None
    record_count = 0
    records = StringIO()

    # Retrieve the auxillary data and extract it
    with open(tape6_filename, "r") as tape6_fd:
        # Skip the beginning of the data that is not needed
        for line in tape6_fd:
            line = line.strip()
            if line.startswith('RADIANCE(WATTS/CM2-STER-XXX)'):
                break

        # This is where our data resides
        for line in tape6_fd:
            # Remove all whitespace and newlines
            line = ' '.join(line.strip().split())

            # Skip warnings, but output them to the log
            if 'WARNING' in line:
                logger.warning("MODTRAN: {0}".format(line))
                continue

            # Skip empty and header lines
            if (not line or
                    line.startswith('RADIANCE') or
                    line.startswith('FREQ') or
                    line.startswith('EMISSION') or
                    line.startswith('(CM-1)')):
                continue

            # Skip the remainding of the file
            if line.startswith('MULTIPLE SCATTERING CALCULATION RESULTS:'):
                break

            parts = line.split(' ')
            if record_count > 0:
                records.write('\n')
            records.write(' '.join([parts[1], parts[12]]))
            record_count += 1

        # Retreive the value from the tape6 file
        target_pixel_surface_temp = extract_tpst(tape6_fd)

    create_output(args, target_pixel_surface_temp, record_count, records)

    records.close()


# ============================================================================
def process_pltout_results(args):
    '''
    Description:
        Parse and clean up the pltout.asc results.
    '''

    logger = logging.getLogger(__name__)

    tape6_filename = '{0}/tape6'.format(args.input_path)
    pltout_filename = '{0}/pltout.asc'.format(args.input_path)

    if not os.path.isfile(tape6_filename):
        raise Exception("Missing tape6 file in input directory")

    if not os.path.isfile(pltout_filename):
        raise Exception("Missing pltout.asc file in input directory")

    target_pixel_surface_temp = None
    record_count = 0
    records = StringIO()

    # Retrieve the auxillary data and extract it
    with open(pltout_filename, "r") as pltout_fd:
        for line in pltout_fd:
            # Remove all whitespace and newlines
            line = ' '.join(line.strip().split())

            if len(line) > 0:
                if record_count > 0:
                    records.write('\n')
                records.write(line)
                record_count += 1

    # Search for the TARGET-PIXEL (H2) SURFACE TEMPERATURE
    with open(tape6_filename, "r") as tape6_fd:
        # Retreive the value from the tape6 file
        target_pixel_surface_temp = extract_tpst(tape6_fd)

    create_output(args, target_pixel_surface_temp, record_count, records)

    records.close()


# ============================================================================
if __name__ == '__main__':
    '''
    Description:
        Performs gathers input parameters and performs the LST MODTRAN results
        processing.
    '''

    # Create a command line arugment parser
    description = ("Reads MODTRAN results from either tape6 or pltout.asc"
                   " files processes them to the follow-on input format"
                   " for the next step in LST generation")
    parser = ArgumentParser(description=description)

    # ---- Add parameters ----
    # Required parameters
    parser.add_argument('--tape6',
                        action='store_true', dest='tape6', required=False,
                        default=False,
                        help="Use the radiance values found in the TAPE6 file")

    parser.add_argument('--pltout',
                        action='store_true', dest='pltout', required=False,
                        default=False,
                        help=("Use the radiance values found in the pltout"
                              " file"))

    parser.add_argument('--input-path',
                        action='store', dest='input_path', required=True,
                        help="Where to find the MODTRAN output files.")

    parser.add_argument('--output-path',
                        action='store', dest='output_path', required=True,
                        help="Where to place the output files.")

    # Parse the command line parameters
    args = parser.parse_args()

    # Setup the default logger format and level. log to STDOUT.
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:'
                                '%(funcName)s -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=logging.INFO,
                        stream=sys.stdout)

    # Get the logger
    logger = logging.getLogger(__name__)

    if ((not args.tape6) and (not args.pltout)):
        logger.fatal("No data source specified.")
        logger.fatal("Error processing LST MODTRAN results."
                     "  Processing will terminate.")
        sys.exit(EXIT_FAILURE)

    if not os.path.isdir(args.input_path):
        logger.fatal("--input-path directory not found")
        sys.exit(EXIT_FAILURE)

    if not os.path.isdir(args.output_path):
        logger.fatal("--output-path directory not found")
        sys.exit(EXIT_FAILURE)

    try:
        if args.tape6:
            logger.info("Using TAPE6 results")
            process_tape6_results(args)
        else:
            logger.info("Using pltout.asc results")
            process_pltout_results(args)

    except Exception, e:
        logger.exception("Error processing LST MODTRAN results."
                         "  Processing will terminate.")
        sys.exit(EXIT_FAILURE)

    logger.info("LST MODTRAN results processing complete")
    sys.exit(EXIT_SUCCESS)
