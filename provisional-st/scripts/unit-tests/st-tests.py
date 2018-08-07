#!/usr/bin/env python


import os
import glob
import unittest
import subprocess
import filecmp
import shutil

class TestST(unittest.TestCase):
    def setUp(self):
        # The XML file for the scene we're working with
        self.xml_filename = 'LE07_L1TP_043028_20020419_20180206_01_T1.xml';

        # List of files linked from ESPA_UNIT_TEST_DATA_DIR into the local
        # unit test directory, to use as inputs for a test
        self.link_files = [
            os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', self.xml_filename),
            ]
        # List of output files and directories to compare to the expected
        # results in ESPA_UNIT_TEST_DATA_DIR
        self.check_files = []
        self.check_dirs = []

    def setUpTestLinks(self):
        for lf in self.link_files:
            if not os.path.exists(os.path.basename(lf)):
                os.symlink(lf, os.path.basename(lf))

    def tearDown(self):
        for lf in self.link_files:
            os.unlink(os.path.basename(lf))
        for d in self.check_dirs:
            shutil.rmtree(os.path.basename(d))
        for cf in self.check_files:
            os.unlink(os.path.basename(cf))

    def test_run(self):
        pass

    def run_test_case(self, cmd):
        status = subprocess.call(cmd)
        self.assertEqual(status, 0)

        for check_dir in self.check_dirs:
            dc = filecmp.dircmp(check_dir, os.path.basename(check_dir))
            self.assertEqual(len(dc.diff_files), 0)
            self.assertEqual(len(dc.funny_files), 0)
            self.assertEqual(len(dc.common_funny), 0)

        for check_file in self.check_files:
            print check_file
            self.assertTrue(
                filecmp.cmp(check_file, os.path.basename(check_file)))

class TestGridPoints(TestST):
    """Test st_determine_grid_points.py """

    def setUp(self):
        super(TestGridPoints, self).setUp()

        self.data_path = os.environ['ST_DATA_DIR']
        self.assertTrue(os.path.exists(self.data_path))

        self.link_files.extend(
            glob.glob(os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', '*bt_band6*')))
        self.setUpTestLinks()

        self.check_files = [
            os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', 'grid_points.bin'),
            os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', 'grid_points.hdr')
            ]

    def test_run(self):
        cmd = ["../st_determine_grid_points.py", "--xml", self.xml_filename,
             "--data_path", self.data_path]
        self.run_test_case(cmd)

class TestAuxMerra(TestST):
    """Test st_extract_auxiliary_merra_data.py """

    def setUp(self):
        super(TestAuxMerra, self).setUp()

        self.aux_path = os.environ['MERRA2_AUX_DIR']
        self.assertTrue(os.path.exists(self.aux_path))

        self.setUpTestLinks()

        self.check_dirs = [
            os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', 'H_t0'),
            os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', 'H_t1'),
            os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', 'QV_t0'),
            os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', 'QV_t1'),
            os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', 'T_t0'),
            os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', 'T_t1')
            ]

    def test_run(self):
        cmd = ["../st_extract_auxiliary_merra_data.py", "--xml", 
               self.xml_filename, "--aux_path", self.aux_path]
        self.run_test_case(cmd)

class TestModtranInput(TestST):
    """Test st_build_modtran_input.py  """

    def setUp(self):
        super(TestModtranInput, self).setUp()

        self.link_files.extend([
            os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', 'grid_points.bin'),
            os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', 'grid_points.hdr'),
            os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', 'H_t0'),
            os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', 'H_t1'),
            os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', 'QV_t0'),
            os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', 'QV_t1'),
            os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', 'T_t0'),
            os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', 'T_t1'),
            ])
        self.link_files.extend(
            glob.glob(os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', 'LE07*elevation*')))
        self.setUpTestLinks()

        self.data_path = os.environ['ST_DATA_DIR']
        self.assertTrue(os.path.exists(self.data_path))

        self.check_dirs = glob.glob(
            os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', '00?_00?_???_???'))

        self.check_files = [
            os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', 'grid_elevations.txt'),
            os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', 'modtran_elevations.txt')
            ]

    def test_run(self):
        cmd = ["../st_build_modtran_input.py", "--xml", self.xml_filename,
             "--data_path", self.data_path]
        self.run_test_case(cmd)

class TestEmissivity(TestST):
    """Test estimate_landsat_emissivity.py and 
            estimate_landsat_emissivity_stdev.py """

    def setUp(self):
        super(TestEmissivity, self).setUp()

        self.link_files.extend(
            glob.glob(os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', '*bt_band6*')))
        self.link_files.extend(
            glob.glob(os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', '*toa_band[2345]*')))
        self.setUpTestLinks()

        self.aster_ged_server = os.environ['ASTER_GED_SERVER']
        self.aster_ged_path = os.environ['ASTER_GED_PATH']
        # Do not test for existence of aster_ged_path, it may be a cloud bucket

        self.check_files = glob.glob(
            os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', 'LE07*_emis*'))
        # Check the tif files from one ASTER granule
        self.check_files.extend(
            glob.glob(os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', 'AG*tif')))
        # Intermediate files
        self.check_files.extend(
            glob.glob(os.path.join(os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
                'espa-surface-temperature', '[a-z]*tif')))

    def test_run(self):
        cmd = ["../estimate_landsat_emissivity.py", "--xml", self.xml_filename,
               "--intermediate",
               "--aster-ged-server-name", self.aster_ged_server,
               "--aster-ged-server-path", self.aster_ged_path]
        self.run_test_case(cmd)

if __name__ == '__main__':
    unittest.main(verbosity=2)
