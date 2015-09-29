'''
    FILE: unit-tests.py

    PURPOSE: Provides unit testing for this directory.

    PROJECT: Land Satellites Data Systems Science Research and Development
             (LSRD) at the USGS EROS

    LICENSE: NASA Open Source Agreement 1.3

    HISTORY:

    Date              Reason
    ----------------  --------------------------------------------------------
    Sep/2015          Initial implementation
'''


import os
import sys
import shutil
import glob
import filecmp
import unittest

# Add the parent directory where the modules to test are located
sys.path.insert(0, '..')
from extract_auxiliary_narr_data import AuxNARRGribProcessor
from lst_environment import Environment


class LSRD_ValidationFramework(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(LSRD_ValidationFramework, self).__init__(*args, **kwargs)
        self.cleanup = True

        if not self.name:
            raise Exception('self.name must be defined')

        # Verify the environment
        self.lsrd_validation_dir = os.environ.get('LSRD_VALIDATION_DIR')
        if self.lsrd_validation_dir is None:
            raise Exception('Missing environment variable LSRD_VALIDATION_DIR')

    def assertFilesEqual(self, file_1, file_2):
        '''Assert that two files are equal or not.'''

        self.cleanup = self.assertTrue(os.path.exists(file_1),
                                       '{0} Does not exist'.format(file_1))
        self.cleanup = self.assertTrue(os.path.exists(file_2),
                                       '{0} Does not exist'.format(file_2))

        self.cleanup = self.assertTrue(filecmp.cmp(file_1, file_2))


class AuxNARRGribProcessor_TestCase(LSRD_ValidationFramework):
    '''Tests for Grib file processing.'''

    def __init__(self, *args, **kwargs):
        self.name = 'AuxNARRGribProcessor_TestCase'
        super(AuxNARRGribProcessor_TestCase, self).__init__(*args, **kwargs)

        # Validation data is presummed to be available if the directory exists
        self.validation_path = os.path.join(self.lsrd_validation_dir,
                                            self.name)
        if not os.path.isdir(self.validation_path):
            raise Exception('Missing validation data for [{0}]'
                            .format(self.name))

        # Define the directories that are produced
        self.directories = ['HGT_1', 'HGT_2',
                            'SPFH_1', 'SPFH_2',
                            'TMP_1', 'TMP_2']

    def setUp(self):
        '''setup'''

        self.input_xml = os.path.join(self.validation_path,
                                      'LT50420342011119PAC01.xml')

        # Specify the XML metadata file defining the data to process
        self.processor = AuxNARRGribProcessor(self.input_xml)

        # Process the associated AUX data
        self.processor.extract_aux_data()

    def tearDown(self):
        '''Cleanup'''

        if self.cleanup:
            for directory in self.directories:

                if os.path.isdir(directory):
                    shutil.rmtree(directory)

    def test_process_grib_data(self):
        '''Test the processing of grib files from our internal archive.'''

        for directory in self.directories:
            self.cleanup = self.assertEqual(True, os.path.isdir(directory))

            # Start with the local files
            files = glob.glob(os.path.join(directory, '*'))

            # Add the validation files
            validation_directory = os.path.join(self.validation_path,
                                                directory)
            files.extend(glob.glob(os.path.join(validation_directory, '*')))

            # We only want the filenames
            files = [os.path.basename(x) for x in files]

            # Make a unique list of the filenames
            files = sorted(list(set(files)))

            # Process through each file
            for filename in files:
                local_file = os.path.join(directory, filename)
                validation_file = os.path.join(validation_directory, filename)

                self.assertFilesEqual(validation_file, local_file)


class Environment_TestCase(LSRD_ValidationFramework):
    '''Tests Environment Class'''

    def __init__(self, *args, **kwargs):
        self.name = 'Environment_TestCase'
        super(Environment_TestCase, self).__init__(*args, **kwargs)

    def setUp(self):
        '''setup'''

        os.environ['LST_DATA_DIR'] = '/usr/local'
        os.environ['LST_AUX_DIR'] = '/usr/local'
        os.environ['ASTER_GED_SERVER_NAME'] = 'ASTER_GED_SERVER_NAME'

        self.environment = Environment()

    def test_LST_DATA_DIR(self):
        '''Test the LST_DATA_DIR environment variable'''

        self.assertEqual('/usr/local',
            self.environment.get_lst_data_directory())

    def test_LST_AUX_DIR(self):
        '''Test the LST_AUX_DIR environment variable'''

        self.assertEqual('/usr/local',
            self.environment.get_lst_aux_directory())

    def test_ASTER_GED_SERVER_NAME(self):
        '''Test the ASTER_GED_SERVER_NAME environment variable'''

        self.assertEqual('ASTER_GED_SERVER_NAME',
            self.environment.get_aster_ged_server_name())


if __name__ == '__main__':
    unittest.main()
