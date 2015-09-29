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


class LSRD_ValidationFramework(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(LSRD_ValidationFramework, self).__init__(*args, **kwargs)
        self.cleanup = True

    def assertFilesEqual(self, file_1, file_2):
        '''Assert that two files are equal or not.'''

        self.cleanup = self.assertTrue(os.path.exists(file_1),
                                       '{0} Does not exist'.format(file_1))
        self.cleanup = self.assertTrue(os.path.exists(file_2),
                                       '{0} Does not exist'.format(file_2))

        self.cleanup = self.assertTrue(filecmp.cmp(file_1, file_2))


class AuxNARRGribProcessor_TestCase(LSRD_ValidationFramework):
    '''Tests for Grib file processing.'''

    def setUp(self):
        '''setup'''

        self.lsrd_validation_dir = os.environ.get('LSRD_VALIDATION_DIR')
        if self.lsrd_validation_dir is None:
            raise Exception('Missing environment variable LSRD_VALIDATION_DIR')

        self.validation_path = os.path.join(self.lsrd_validation_dir,
                                            'AuxNARRGribProcessor_TestCase')
        self.input_xml = os.path.join(self.validation_path,
                                      'LT50420342011119PAC01.xml')

        # Specify the XML metadata file defining the data to process
        self.processor = AuxNARRGribProcessor(self.input_xml)

        # Define the directories that will be produced
        self.directories = ['HGT_1', 'HGT_2',
                            'SPFH_1', 'SPFH_2',
                            'TMP_1', 'TMP_2']

        # Define the files that will be produced in each directory
        self.files = ['100.txt', '125.txt', '150.txt', '175.txt',
                      '200.txt', '225.txt', '250.txt', '275.txt',
                      '300.txt', '350.txt',
                      '400.txt', '450.txt',
                      '500.txt', '550.txt',
                      '600.txt', '650.txt',
                      '700.txt', '725.txt', '750.txt', '775.txt',
                      '800.txt', '825.txt', '850.txt', '875.txt',
                      '900.txt', '925.txt', '950.txt', '975.txt',
                      '1000.txt']

    def tearDown(self):
        '''Cleanup'''

        if self.cleanup:
            for directory in self.directories:

                if os.path.isdir(directory):
                    shutil.rmtree(directory)

    def test_process_grib_data(self):
        '''Test the processing of grib files from our internal archive.'''

        self.processor.extract_aux_data()

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


class Next_TestCase(LSRD_ValidationFramework):
    '''Tests for XXXXX file processing.'''

    def setUp(self):
        '''setup'''

        pass

    def tearDown(self):
        '''Cleanup'''

        if self.cleanup:
            # Add your cleanup code here
            pass

    def test_something(self):
        '''Something'''

        pass


if __name__ == '__main__':
    unittest.main()
