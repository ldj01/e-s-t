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
import unittest

# Add the parent directory where the modules to test are located
sys.path.insert(0, '..')
from extract_auxiliary_narr_data import AuxNARRGribProcessor


class AuxNARRGribProcessor_TestCase(unittest.TestCase):
    '''Tests for Grib file processing.'''

    def setUp(self):
        '''setup'''

        # Specify the XML metadata file defining the data to process
        self.processor = AuxNARRGribProcessor('LT50420342011119PAC01.xml')

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

        for directory in self.directories:

            for filename in self.files:
                path = os.path.join(directory, filename)
                if os.path.exists(path):
                    os.unlink(path)

            if os.path.isdir(directory):
                os.rmdir(directory)

    def test_process_grib_data(self):
        '''Test the processing of grib files from our internal archive.'''

        self.processor.extract_aux_data()

        for directory in self.directories:
            self.assertEqual(True, os.path.isdir(directory))

            for filename in self.files:
                path = os.path.join(directory, filename)
                self.assertEqual(True, os.path.exists(path))


if __name__ == '__main__':
    unittest.main()
