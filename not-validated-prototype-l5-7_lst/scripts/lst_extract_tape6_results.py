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
def process_tape6_results(args):
    '''
    Description:
        Parses the tape6 and extracts the WAVELENGTH and TOTAL_RADIANCE
        values.  Writes the values to the requested output file.
    '''

    logger = logging.getLogger(__name__)

    results_buffer = StringIO()

    # Retrieve the auxillary data and extract it
    with open(args.tape6_filename, "r") as tape6_fd:
        # Skip the beginning of the data that is not needed
        for line in tape6_fd:
            line = line.strip()
            if line.startswith('RADIANCE(WATTS/CM2-STER-XXX)'):
                break

        for line in tape6_fd:
            # Remove all whitespace and newlines
            line = ' '.join(line.strip().split())

            # Skip warnings, but output them to the log
            if 'WARNING' in line:
                logger.warning("MODTRAN: {0}".format(line))
                continue

            # Skip empty and header lines
            if (not line
                    or line.startswith('RADIANCE')
                    or line.startswith('FREQ')
                    or line.startswith('EMISSION')
                    or line.startswith('(CM-1)')):
                continue

            # Skip the remainding of the file
            if line.startswith('MULTIPLE SCATTERING CALCULATION RESULTS:'):
                break

            parts = line.split(' ')
            results_buffer.write(' '.join([parts[1], parts[12]]))
            results_buffer.write('\n')

    with open(args.parsed_filename, "w") as parsed_fd:
        parsed_fd.write(results_buffer.getvalue())

    results_buffer.close()


# ============================================================================
def process_pltout_results(args):
    '''
    Description:
        Parse and clean up the pltout.asc results.
    '''

    logger = logging.getLogger(__name__)

    results_buffer = StringIO()

    # Retrieve the auxillary data and extract it
    with open(args.pltout_filename, "r") as pltout_fd:
        for line in pltout_fd:
            # Remove all whitespace and newlines
            line = ' '.join(line.strip().split())

            results_buffer.write(line)
            results_buffer.write('\n')

    with open(args.parsed_filename, "w") as parsed_fd:
        parsed_fd.write(results_buffer.getvalue())

    results_buffer.close()


# ============================================================================
if __name__ == '__main__':
    '''
    Description:
        Performs gathers input parameters and performs the LST tape6
        processing.
    '''

    # Create a command line arugment parser
    description = ("Retrieves and generates auxillary LST inputs, then"
                   " processes and calls other executables for LST generation")
    parser = ArgumentParser(description=description)

    # ---- Add parameters ----
    # Required parameters
    parser.add_argument('--tape6',
                        action='store', dest='tape6_filename', required=False,
                        default=None,
                        help="The TAPE6 file to process")

    parser.add_argument('--pltout',
                        action='store', dest='pltout_filename', required=False,
                        default=None,
                        help="The output parsed results file")

    parser.add_argument('--parsed',
                        action='store', dest='parsed_filename', required=True,
                        help="The output parsed results file")

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

    try:
        if not args.pltout_filename:
            if args.tape6_filename == '':
                logger.fatal("No TAPE6 filename provided.")
                logger.fatal("Error processing LST TAPE6 results."
                             "  Processing will terminate.")
                sys.exit(EXIT_FAILURE)

            logger.info("Using TAPE6 results")
            process_tape6_results(args)
        else:
            if args.pltout_filename == '':
                logger.fatal("No pltout filename provided.")
                logger.fatal("Error processing LST pltout results."
                             "  Processing will terminate.")
                sys.exit(EXIT_FAILURE)

            logger.info("Using pltout.asc results")
            process_pltout_results(args)

    except Exception, e:
        logger.exception("Error processing LST TAPE6 results."
                         "  Processing will terminate.")
        sys.exit(EXIT_FAILURE)

    logger.info("LST TAPE6 results processing complete")
    sys.exit(EXIT_SUCCESS)
