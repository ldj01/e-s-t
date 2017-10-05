#! /usr/bin/env python

'''
    PURPOSE: Retrieves archived NARR files from the NCEP for the dates
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
import shutil
import logging
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from datetime import datetime, timedelta, date
import collections
import st_aux_config as config

from st_aux_utilities import System
from st_aux_http_session import HttpSession
from st_aux_version import VERSION_TEXT
from st_aux_parameters import NARR_VARIABLES


class Ncep(object):
    '''Interface for interacting with NCEP website

    NarrData depends on the following functionality of this class:
        Provide means of getting a particular grib file into current working
            directory via get_grib_file()
        Provides a dict of all files available via get_dict_of_date_modified()
            keys are filenames, values are time of last modification
        Provides format of the grib filename via get_filename()
    '''
    mtime_by_name = None
    session = None

    @staticmethod
    def get_url(cfg, filename):
        '''Return the URL for external retrieval of the file'''
        return cfg.ncep.data_url_format.format(filename)

    @staticmethod
    def get_filename(cfg, year, month, day, hour):
        '''Return the filename to grab on the external system'''
        fmt = cfg.ncep.data_name_format
        return fmt.format(year, month, day, hour)

    @staticmethod
    def get_datetime_from_filename(filename):
        '''Extracts tuple (year, month, day, hour) from filename

        Precondition:
            File is of format "rcdas.2015081321.awip32.merged" for the
                measurement at 21:00 (9pm) on August 13th and similar for all
                other measurements
        Postcondition:
            Returns datetime object of the time when the measurement took place
        '''
        return datetime.strptime(filename.split('.')[1],
                                 '%Y%m%d%H')

    @classmethod
    def get_grib_file(cls, cfg, filename):
        '''Retrieves the grib file from the external system

        Precondition:
            File with "filename" exists on the NCEP website.
            "get_url" will return address of file on website
            "get_external_last_modified" returns a datetime object
                This datetime object should be the last time the item was
                updated on the website.
            File does not exists in current directory
        Postcondition:
            grib_file associated with this filename is in current directory
            Logs the last modified time of file and the address used.
            If file already exists then only an info message will be recorded.
        '''
        logger = logging.getLogger(__name__)

        if os.path.isfile(filename):
            logger.info('{0} already exists. Skipping download.'
                        .format(filename))
        else:
            logger.info('Retrieving {0}'.format(filename))
            cls.get_session(cfg).http_transfer_file(cls.get_url(cfg, filename),
                                                 filename)

    @classmethod
    def get_list_of_external_data(cls, cfg):
        '''Retrieves list of available data from website's directory listing

        Sample line from url reqest the list of files (single line):
        '<tr><td><a href="rcdas.2015010300.awip32.merged.b">
            rcdas.2015010300.awip32.merged.b</a></td>
        <td align="right">08-Jan-2015 10:12  </td>
        <td align="right">1.3M</td></tr>\n'
        '''
        archive_data = collections.namedtuple('ArchiveData',
                                              ['name', 'mtime', 'size'])

        lines_thrown = 0
        data_list = list()

        custom_session = cls.get_session(cfg)
        try:
            data = custom_session.get_lines_from_url(Ncep.get_url(cfg, ''))

            for line in data:
                if 'awip' not in line:
                    lines_thrown = lines_thrown + 1
                    continue  # go to next line

                line_parts = line.split('">')
                name = line_parts[0].split('="')[1]
                mtime = line_parts[2].split('<')[0]
                size = line_parts[3].split('<')[0]

                mtime = mtime.strip()  # Remove extra space
                size = size.strip()  # Remove extra space
                data_list.append(archive_data(name=name,
                                              mtime=mtime,
                                              size=size))

        except Exception:
            raise

        return data_list

    @classmethod
    def get_session(cls, cfg):
        '''Obtains and then retains session used for downloading'''

        if cls.session is None:
            # Establish an HTTP session
            cls.session = HttpSession()

        return cls.session

    @classmethod
    def get_dict_of_date_modified(cls, cfg):
        '''Returns a dictionary of mtime for ext. files with filename as key

        Note:
            If the dictionary has been cached then use it otherwise create it.
        Precondition:
            Requires that get_list_of_external_data() returns list containing
                NamedTuples with tuple.mtime and tuple.name defined.
        Postcondition: Returns a dictionary
            filename as key, external last modified time as value
        '''

        if cls.mtime_by_name is None:
            data_list = Ncep.get_list_of_external_data(cfg)
            cls.mtime_by_name = {}
            for item in data_list:
                date_modified = datetime.strptime(item.mtime,
                                                  '%d-%b-%Y %H:%M')
                cls.mtime_by_name[item.name] = date_modified

        return cls.mtime_by_name


class NarrData(object):
    '''Interface for interacting with NARR data'''

    class FileMissing(Exception):
        '''Exception raised when file is missing internally or on website'''
        pass

    def __init__(self, year, month, day, hour=00):
        hour = hour/3*3  # Ensures it is a multiple of 3
        self.dt = datetime(year, month, day, hour=hour)

    @staticmethod
    def get_next_narr_data_gen(s_date, e_date, interval=timedelta(hours=3)):
        '''Generator to iterate through NarrData objects in time interval

        Note: time of datetime objects is ignored.
            start of day is used for s_date and end of day for e_date.
        Precondition:
            s_date and e_date are of type datetime.datetime
            s_date and e_date can also be of type datetime.date
            interval is of type datetime.timedelta
        Postcondition: Returns NarrData instance
            s_date < datetime of NarrData < e_date
            Consider narr[i-1], narr[i], narr[i+1] as previous, current and
                next values of iterator. The following statement will be true:
                narr[i-1].dt + interval == narr[i].dt == narr[i].dt - interval
        '''

        try:  # Handles if datetime objects are passed in.
            start_time = datetime.combine(s_date.date(), datetime.min.time())
            end_time = datetime.combine(e_date.date(), datetime.max.time())
        except AttributeError:  # Handles if date objects are passed in.
            start_time = datetime.combine(s_date, datetime.min.time())
            end_time = datetime.combine(e_date, datetime.max.time())

        logger = logging.getLogger(__name__)
        logger.info('Generating list of data from {0} to {1}'
                    .format(start_time.isoformat(), end_time.isoformat()))

        current = NarrData(year=start_time.year, month=start_time.month,
                           day=start_time.day, hour=start_time.hour)
        while current.dt <= end_time:
            yield current
            current = current.get_next(interval)

    def get_internal_directory(self, cfg):
        '''Returns the internal path to the archive'''
        return NarrArchive.get_arch_dir(cfg, self.dt.year, self.dt.month,
                                        self.dt.day)

    def get_internal_filename(self, cfg, variable, ext):
        '''Returns an internally formatted archive filename'''
        return NarrArchive.get_arch_filename(cfg, variable, self.dt.year,
                                             self.dt.month, self.dt.day,
                                             self.dt.hour, ext)

    def get_internal_last_modified(self, cfg, variable='HGT', ext='hdr'):
        '''Stat internal file for mtime. Default to HGT's hdr file.

        Precondition:
            File must exist at path given by get_internal_filepath
        Postcondition:
            returns time of last modification of internal file
            raises NarrData.FileMissing if precondition is violated
        '''
        try:
            filepath = os.path.join(self.get_internal_directory(cfg),
                                    self.get_internal_filename(cfg, variable, 
                                                               ext))
            ts_epoch = os.stat(filepath).st_mtime
            mtime = datetime.fromtimestamp(ts_epoch)
        except OSError:  # Expecting 'No such file or directory'
            raise NarrData.FileMissing

        return mtime

    def get_external_filename(self, cfg):
        '''Returns the name of the grib file as choosen by data source'''
        return Ncep.get_filename(cfg, self.dt.year, self.dt.month, self.dt.day,
                                 self.dt.hour)

    def get_external_last_modified(self, cfg):
        '''Returns last_modified time from dictionary stored in Ncep

        Precondition:
            Ncep class must be able to obtain list of data from website.
            filename must exist as key in dict of date modified
        Postcondition:
            returns date of last modification of the entry with filename as key
            Raises NarrData.FileMissing if either precondition is violated.
        '''
        filename = self.get_external_filename(cfg)

        try:
            # Last modified time according to http table on website
            table_last_mod = Ncep.get_dict_of_date_modified(cfg)[filename]
        except KeyError:
            raise NarrData.FileMissing

        return table_last_mod

    def need_to_update(self, cfg):
        '''Returns boolean of whether file neads to be downloaded

        Precondition:
            get_internal_last_modified and get_external_last_modified return
                value of type datetime.
        Postcondition:
            returns True if either of these conditions are true:
                (1) Internal copy of file is older than remote copy
                (2) Internal file is missing
            returns False if either of these conditions are true:
                (1) Internal modified time is more recent then external
                    modified time
                (2) External file is missing
        '''
        logger = logging.getLogger(__name__)
        try:
            ext_mtime = self.get_external_last_modified(cfg)
        except NarrData.FileMissing:  # Expecting 'No such file or directory'
            logger.debug('{0} is missing from list of external files'
                         .format(self.get_external_filename(cfg)))
            return False  # File is not available to download

        try:
            # Check if existing data is stale
            return self.get_internal_last_modified(cfg) < ext_mtime
        except NarrData.FileMissing:  # Expecting 'No such file or directory'
            return True  # The file does not exist internally.

    def get_grib_file(self, cfg):
        '''retrieves the grib file'''
        Ncep.get_grib_file(cfg, self.get_external_filename(cfg))

    def extract_vars_from_grib(self, cfg):
        '''process_grib_for_variable for each var in NARR_VARIABLES'''
        for var in NARR_VARIABLES:
            self.process_grib_for_variable(cfg, var)

    def move_files_to_archive(self, cfg):
        '''move_to_archive for each var in NARR_VARIABLES'''
        for var in NARR_VARIABLES:
            self.move_to_archive(cfg, var)

    def remove_grib_file(self, cfg):
        '''removes the grib file'''
        logger = logging.getLogger(__name__)
        logger.debug('ExternalFile(Exists:{0}, Name:{1})'
                     .format(os.path.exists(self.get_external_filename(cfg)),
                             self.get_external_filename(cfg)))
        if os.path.exists(self.get_external_filename(cfg)):
            os.unlink(self.get_external_filename(cfg))

    def get_next(self, time_increment=timedelta(hours=3)):
        '''returns the next NarrData object'''
        next_date = self.dt + time_increment
        return NarrData(year=next_date.year, month=next_date.month,
                        day=next_date.day, hour=next_date.hour)

    def process_grib_for_variable(self, cfg, variable, verbose=False):
        '''Extract the specified variable from the grib file and archive it.

        Precondition:
            A grib file, with the name get_external_filename(), exists in
                current working directory.
            wgrib must be installed on the system
        Postcondition:
            A grib and header for variable will exist in current working
                directory with the name given by get_internal_filename()
        '''
        logger = logging.getLogger(__name__)

        grib_file = self.get_external_filename(cfg)
        hdr_name = self.get_internal_filename(cfg, variable, 'hdr')
        grb_name = self.get_internal_filename(cfg, variable, 'grb')

        if os.path.isfile(grb_name) and os.path.isfile(hdr_name):
            logger.info('{0} and {1} already exist. Skipping extraction.'
                        .format(hdr_name, grb_name))
            return
        logger.info("Processing [{0}]".format(grib_file))

        # Create inventory/header file to extract the variable data
        cmd_create_temp_header = ['wgrib', grib_file, '|', 'grep', variable,
                                  '>', hdr_name]
        # Create grib files for each variable
        cmd_create_final_grib = ['cat', hdr_name, '|',
                                 'wgrib', grib_file, '-i', '-grib',
                                 '-o', grb_name]
        # Create new inventory/header file for the variable
        cmd_create_final_header = ['wgrib', grb_name, '|', 'grep', variable,
                                   '>', hdr_name]

        cmds = [cmd_create_temp_header, cmd_create_final_grib,
                cmd_create_final_header]

        for cmd_list in cmds:
            cmd = ' '.join(cmd_list)
            output = ''
            logger.info('Executing [{0}]'.format(cmd))
            output = System.execute_cmd(cmd)
            if verbose:
                if len(output) > 0:
                    logger.info(output)

    def move_to_archive(self, cfg, variable):
        '''Moves grb and hdr files to archive location.

        Precondition:
            Header and Grib files for variable exist in current directory
        Postcondition:
            Header and Grib files for variable exist in archive directory
            Header and Grib files for variable don't exist in current directory
        '''
        logger = logging.getLogger(__name__)

        dest_path = self.get_internal_directory(cfg)  # Determine the directory
        hdr_name = self.get_internal_filename(cfg, variable, 'hdr')
        grb_name = self.get_internal_filename(cfg, variable, 'grb')

        System.create_directory(dest_path)  # create it if it does not exist

        # Archive the files
        logger.info('Archiving into [{0}]'.format(dest_path))
        # GRIB
        dest_file = os.path.join(dest_path, grb_name)
        shutil.copyfile(grb_name, dest_file)
        # HEADER
        dest_file = os.path.join(dest_path, hdr_name)
        shutil.copyfile(hdr_name, dest_file)

        # Cleanup the working directory
        if os.path.exists(grb_name):
            os.unlink(grb_name)
        if os.path.exists(hdr_name):
            os.unlink(hdr_name)


class NarrArchive(object):
    '''Class for storing the NARR archive data'''
    _base_aux_dir = None

    @classmethod
    def get_arch_filename(cls, cfg, variable, year, month, day, hour, ext):
        '''Determine the filename for the archive data'''
        return (cfg.archive_name_format
                .format(variable, year, month, day, hour*100, ext))

    @classmethod
    def get_arch_dir(cls, cfg, year, month, day):
        '''Determine the directory name for the archive data'''
        return (cfg.archive_directory_format
                .format(cfg.base_archive_directory, year, month, day))


def update(cfg, data_to_be_updated):
    '''Downloads, extracts vars, and cleans temp files for data passed in

    Precondition:
        data_to_be_updated is a list of NarrData objects
        External files exist for every data item
    Postcondition:
        Header and Grib files for each variable for each item in
            data_to_be_updated exist in archive directory
        No temporary files exist in the working directory
    '''

    for data in data_to_be_updated:
        try:
            data.get_grib_file(cfg)
            data.extract_vars_from_grib(cfg)
            data.move_files_to_archive(cfg)
        finally:
            data.remove_grib_file(cfg)


def report(cfg, data_to_report):
    '''Provides measured time, internal mtime and external mtime of data
       passed in

    Note:
        Reports number of files to be downloaded
        Includes header to describe data being output.
        Reports [measured time, internal mtime, external mtime] as csv
    '''

    report_msg = []
    report_msg.append('Measured, Local TimeStamp, Remote TimeStamp')  # Header

    for data in data_to_report:
        line = []
        line.append(data.dt.isoformat())  # Measured datetime

        try:
            line.append(data.get_internal_last_modified(cfg).isoformat())
        except NarrData.FileMissing:
            line.append('-')

        try:
            line.append(data.get_external_last_modified(cfg).isoformat())
        except NarrData.FileMissing:
            line.append('-')

        report_msg.append(', '.join(line))

    print ('\n'.join(report_msg))


def arg_date_type(datestring):
    '''Validates the date string to be a specified format'''

    try:
        return datetime.strptime(datestring, '%Y%m%d').date()
    except ValueError:
        print('Dates must be the in the format: "YYYYMMDD"')
        raise


def parse_arguments(cfg):
    '''
    Description:
        Parses arguments from the command line.
    '''

    default_date_range = int(cfg.search_date_range)

    # Create a command line arugment parser
    description = ('Downloads ST auxiliary inputs, then archives them for'
                   ' future use. Dates must be the in the format: "YYYYMMDD"')
    parser = ArgumentParser(description=description,
                            formatter_class=ArgumentDefaultsHelpFormatter)

    # ---- Add parameters ----
    parser.add_argument('--version',
                        action='version',
                        version=VERSION_TEXT)

    parser.add_argument('--start-date',
                        action='store', dest='start_date',
                        metavar='YYYYMMDD', type=arg_date_type,
                        required=False,
                        default=(date.today() -
                                 timedelta(days=default_date_range)),
                        help='The start date of the date range of auxiliary'
                             ' data to download.')

    parser.add_argument('--end-date',
                        action='store', dest='end_date',
                        metavar='YYYYMMDD', type=arg_date_type,
                        required=False,
                        default=date.today(),
                        help='The end date of the date range of auxiliary'
                             ' data to download.')

    parser.add_argument('--date',
                        action='store', dest='date',
                        metavar='YYYYMMDD', type=arg_date_type,
                        required=False,
                        help='Sets both start and end date to this date.'
                             ' Overrides start-date and end-date arguments.')

    parser.add_argument('--report',
                        action='store_true', dest='report',
                        default=False,
                        help='Only report what will happen.')

    parser.add_argument('--verbose',
                        action='store_true', dest='verbose',
                        default=False,
                        help='Turn verbose logging on.')

    parser.add_argument('--debug',
                        action='store_true', dest='debug',
                        default=False,
                        help='Turn debug logging on.')

    # Parse the command line parameters
    args = parser.parse_args()

    # Check if date was specified. If so then override start and end.
    if args.date is not None:
        args.start_date = args.date
        args.end_date = args.date

    return args


def setup_logging(debug_logging, verbose_logging):
    '''
    Description:
        Configures the logging for this script based on some command line
        parameters.
    '''

    log_level = logging.WARN
    if debug_logging:
        log_level = logging.DEBUG
    elif verbose_logging:
        log_level = logging.INFO

    # Setup the default logger format and level. log to STDOUT
    logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                ' %(levelname)-8s'
                                ' %(filename)s:%(lineno)d:'
                                '%(funcName)s -- %(message)s'),
                        datefmt='%Y-%m-%d %H:%M:%S',
                        level=log_level,
                        stream=sys.stdout)

    # Turn down the requests and urllib3 logging
    logging.getLogger("requests").setLevel(logging.WARNING)
    logging.getLogger("urllib3").setLevel(logging.WARNING)


def main():
    '''
    Description:
        Ensures all data between start_date and end_date are up to date.

    Precondition:
        start_date and end_date are of type datetime.datetime
        start_date and end_date can also be of type datetime.date
    '''

    # Read the configuration file
    cfg = config.get_config()

    # Parse the command-line arguments 
    cmd_args = parse_arguments(cfg)

    # Setup logging 
    setup_logging(cmd_args.debug, cmd_args.verbose)

    logger = logging.getLogger(__name__)

    try:
        # Determine the data that exists within the date range
        data = NarrData.get_next_narr_data_gen(cmd_args.start_date,
                                               cmd_args.end_date)

        # Determine which files are stale or missing internally.
        data_to_be_updated = filter(lambda x: x.need_to_update(cfg), data)
        if len(data_to_be_updated) == 0:
            logger.info('No data found for updating archive')
        else:
            logger.info('Will download {0} files'.
                        format(len(data_to_be_updated)))
        if cmd_args.report:
            report(cfg, list(data_to_be_updated))
        else:
            update(cfg, data_to_be_updated)

    except Exception:
        logger.exception('Processing Failed')
        sys.exit(1)  # EXIT FAILURE

    sys.exit(0)  # EXIT SUCCESS

if __name__ == '__main__':
    main()
