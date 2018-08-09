#!/usr/bin/env python


import os
import glob
import unittest
import subprocess
import filecmp
import shutil
import sys
import logging

# Import the ST modules that are not standalone scripts
sys.path.insert(0, '..')
import build_st_data

class TestST(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(TestST, self).__init__(*args, **kwargs)
        self.unit_test_data_dir = os.path.join(
            os.environ['ESPA_UNIT_TEST_DATA_DIR'], 
            'espa-surface-temperature')

    def setUp(self):
        # The XML file for the scene we're working with
        self.xml_filename = 'LE07_L1TP_043028_20020419_20180206_01_T1.xml';

        # List of files linked from ESPA_UNIT_TEST_DATA_DIR into the local
        # unit test directory, to use as inputs for a test
        self.link_files = [
            os.path.join(self.unit_test_data_dir, self.xml_filename),
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
        self.do_checks()

    def do_checks(self):
        for check_dir in self.check_dirs:
            dc = filecmp.dircmp(check_dir, os.path.basename(check_dir))
            self.assertEqual(len(dc.diff_files), 0)
            self.assertEqual(len(dc.funny_files), 0)
            self.assertEqual(len(dc.common_funny), 0)

        for check_file in self.check_files:
            print "Comparing", check_file
            self.assertTrue(
                filecmp.cmp(check_file, os.path.basename(check_file)))

    def compare_updated_xml(self, extension):
        # Compare the updated XML file with the expected one, which is 
        # named with the given extension
        expected_xml_filename = os.path.join(self.unit_test_data_dir, 
                self.xml_filename + extension)
        print "Comparing", self.xml_filename, "to", expected_xml_filename
        # Ignore app_version and production_date
        diff_cmd = ["diff", "-I", "production_date", "-I", "app_version", 
            self.xml_filename, expected_xml_filename]
        status = subprocess.call(diff_cmd)
        self.assertEqual(status, 0)
        os.unlink(self.xml_filename)

class TestGridPoints(TestST):
    """Test st_determine_grid_points.py """

    def setUp(self):
        super(TestGridPoints, self).setUp()

        self.data_path = os.environ['ST_DATA_DIR']
        self.assertTrue(os.path.exists(self.data_path))

        self.link_files.extend(
            glob.glob(os.path.join(self.unit_test_data_dir, '*bt_band6*')))
        self.setUpTestLinks()

        self.check_files = [
            os.path.join(self.unit_test_data_dir, 'grid_points.bin'),
            os.path.join(self.unit_test_data_dir, 'grid_points.hdr')
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
            os.path.join(self.unit_test_data_dir, 'H_t0'),
            os.path.join(self.unit_test_data_dir, 'H_t1'),
            os.path.join(self.unit_test_data_dir, 'QV_t0'),
            os.path.join(self.unit_test_data_dir, 'QV_t1'),
            os.path.join(self.unit_test_data_dir, 'T_t0'),
            os.path.join(self.unit_test_data_dir, 'T_t1')
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
            os.path.join(self.unit_test_data_dir, 'grid_points.bin'),
            os.path.join(self.unit_test_data_dir, 'grid_points.hdr'),
            os.path.join(self.unit_test_data_dir, 'H_t0'),
            os.path.join(self.unit_test_data_dir, 'H_t1'),
            os.path.join(self.unit_test_data_dir, 'QV_t0'),
            os.path.join(self.unit_test_data_dir, 'QV_t1'),
            os.path.join(self.unit_test_data_dir, 'T_t0'),
            os.path.join(self.unit_test_data_dir, 'T_t1'),
            ])
        self.link_files.extend(
            glob.glob(os.path.join(self.unit_test_data_dir, 'LE07*elevation*')))
        self.setUpTestLinks()

        self.data_path = os.environ['ST_DATA_DIR']
        self.assertTrue(os.path.exists(self.data_path))

        self.check_dirs = glob.glob(
            os.path.join(self.unit_test_data_dir, '00?_00?_???_???'))

        self.check_files = [
            os.path.join(self.unit_test_data_dir, 'grid_elevations.txt'),
            os.path.join(self.unit_test_data_dir, 'modtran_elevations.txt')
            ]

    def test_run(self):
        cmd = ["../st_build_modtran_input.py", "--xml", self.xml_filename,
             "--data_path", self.data_path]
        self.run_test_case(cmd)

class TestEmissivity(TestST):
    """Test estimate_landsat_emissivity.py """

    def setUp(self):
        super(TestEmissivity, self).setUp()

        # This test modifies the XML file, so instead of linking to it, 
        # make a copy.  Also remove it from the link_files list.
        shutil.copyfile(
            os.path.join(self.unit_test_data_dir, self.xml_filename), 
            self.xml_filename)
        self.link_files = []

        self.link_files.extend(
            glob.glob(os.path.join(self.unit_test_data_dir, '*bt_band6*')))
        self.link_files.extend(
            glob.glob(os.path.join(self.unit_test_data_dir, 
                '*toa_band[2345]*')))
        self.setUpTestLinks()

        self.aster_ged_server = os.environ['ASTER_GED_SERVER']
        self.aster_ged_path = os.environ['ASTER_GED_PATH']
        # Do not test for existence of aster_ged_path, it may be a cloud bucket

        self.check_files = glob.glob(
            os.path.join(self.unit_test_data_dir, 'LE07*_emis.*'))
        # Check the tif files from one ASTER granule
        self.check_files.extend(
            glob.glob(os.path.join(self.unit_test_data_dir, 'AG*emis.tif')))
        self.check_files.extend(
            glob.glob(os.path.join(self.unit_test_data_dir, 'AG*ndvi.tif')))
        # Intermediate files
        self.check_files.extend(
            glob.glob(os.path.join(self.unit_test_data_dir, '[a-z]*ndvi*.tif')))
        self.check_files.extend(
            glob.glob(os.path.join(self.unit_test_data_dir, 
                'landsat_emis_[mw]*.tif')))

    def test_run(self):
        cmd = ["../estimate_landsat_emissivity.py", "--xml", self.xml_filename,
               "--intermediate",
               "--aster-ged-server-name", self.aster_ged_server,
               "--aster-ged-server-path", self.aster_ged_path]
        self.run_test_case(cmd)
        self.compare_updated_xml(".emis")
        # Clean up any leftover ASTER files, not already in the check_files list
        check_files_base = map(lambda x:os.path.basename(x), self.check_files)
        for aster_file in glob.glob("AG100.v003*"):
            if check_files_base.count(aster_file) > 0:
                continue
            os.unlink(aster_file)

class TestEmissivityStdev(TestST):
    """Test estimate_landsat_emissivity_stdev.py """

    def setUp(self):
        super(TestEmissivityStdev, self).setUp()

        # This test modifies the XML file, so instead of linking to it, 
        # make a copy.  Also remove it from the link_files list.
        shutil.copyfile(
            os.path.join(self.unit_test_data_dir, self.xml_filename), 
            self.xml_filename)
        self.link_files = []

        self.link_files.extend(
            glob.glob(os.path.join(self.unit_test_data_dir, '*bt_band6*')))
        self.link_files.extend(
            glob.glob(os.path.join(self.unit_test_data_dir, '*toa_band3*')))
        self.setUpTestLinks()

        self.aster_ged_server = os.environ['ASTER_GED_SERVER']
        self.aster_ged_path = os.environ['ASTER_GED_PATH']
        # Do not test for existence of aster_ged_path, it may be a cloud bucket

        self.check_files = glob.glob(
            os.path.join(self.unit_test_data_dir, 'LE07*_emis_*'))
        # Check the tif files from one ASTER granule
        self.check_files.extend(
            glob.glob(os.path.join(self.unit_test_data_dir, 
                'AG*emis_stdev.tif')))
        # Intermediate files
        self.check_files.extend(
            glob.glob(os.path.join(self.unit_test_data_dir, 
                'landsat_emis_stdev*.tif')))

    def test_run(self):
        cmd = ["../estimate_landsat_emissivity_stdev.py", "--xml", 
                self.xml_filename, "--intermediate",
               "--aster-ged-server-name", self.aster_ged_server,
               "--aster-ged-server-path", self.aster_ged_path]
        self.run_test_case(cmd)
        self.compare_updated_xml(".emis_stdev")
        # Clean up any leftover ASTER files, not already in the check_files list
        check_files_base = map(lambda x:os.path.basename(x), self.check_files)
        for aster_file in glob.glob("AG100.v003*"):
            if check_files_base.count(aster_file) > 0:
                continue
            os.unlink(aster_file)

class TestModtran(TestST):
    """Test st_run_modtran.py  """

    def setUp(self):
        super(TestModtran, self).setUp()

        self.link_files.extend([
            os.path.join(self.unit_test_data_dir, 'grid_points.bin'),
            os.path.join(self.unit_test_data_dir, 'grid_points.hdr'),
            ])
        self.setUpTestLinks()

        # Modtran output files are created in the point directories
        # So set up special list of output files to check in the subdir
        self.check_subdir_files = glob.glob(
            os.path.join(self.unit_test_data_dir, '*', '*', '*', '*', 
                'st_modtran.data'))
        self.check_subdir_files.extend(
            glob.glob(os.path.join(self.unit_test_data_dir, '*', '*', '*', '*', 
                'st_modtran.hdr')))

        # Modtran runs and creates output directly in the point directories
        # So set up copies of them
        self.copy_dirs = glob.glob(
            os.path.join(self.unit_test_data_dir, '00?_00?_???_???'))
        for d in self.copy_dirs:
            shutil.copytree(d, os.path.basename(d))

        self.modtran_data_path = os.environ['MODTRAN_DATA_DIR']
        self.assertTrue(os.path.exists(self.modtran_data_path))

        if os.environ.has_key('NUM_THREADS'):
            self.num_threads = os.environ['NUM_THREADS']
        else:
            self.num_threads = 4

    def tearDown(self):
        super(TestModtran, self).tearDown()
        # Delete extra directory copies
        for d in self.copy_dirs:
            shutil.rmtree(os.path.basename(d))

    def test_run(self):
        cmd = ["../st_run_modtran.py",
               "--process_count", str(self.num_threads),
               "--modtran_data_path", self.modtran_data_path]
        self.run_test_case(cmd)

        # Compare the output that was created in the point directories
        for check_subdir_file in self.check_subdir_files:
            local_subdir_file = check_subdir_file.replace(
                os.path.join(self.unit_test_data_dir, ""))
            print "Comparing", local_subdir_file
            self.assertTrue(
                filecmp.cmp(check_subdir_file, local_subdir_file));

class TestAtmosParam(TestST):
    """Test st_atmospheric_parameters (C code) """

    def setUp(self):
        super(TestAtmosParam, self).setUp()

        self.data_path = os.environ['ST_DATA_DIR']
        self.assertTrue(os.path.exists(self.data_path))

        # This test modifies the XML file, so instead of linking to it, 
        # make a copy.  Also remove it from the link_files list.
        shutil.copyfile(
            os.path.join(self.unit_test_data_dir, self.xml_filename), 
            self.xml_filename)
        self.link_files = []

        self.link_files.extend(
            glob.glob(os.path.join(self.unit_test_data_dir, 'LE07*b61*')))
        self.link_files.extend(
            glob.glob(os.path.join(self.unit_test_data_dir, 'LE07*elevation*')))
        self.link_files.extend(
            glob.glob(os.path.join(self.unit_test_data_dir, '00?_00?_???_???')))
        self.link_files.extend([
            os.path.join(self.unit_test_data_dir, 'grid_points.bin'),
            os.path.join(self.unit_test_data_dir, 'grid_points.hdr'),
            os.path.join(self.unit_test_data_dir, 'grid_elevations.txt'),
            os.path.join(self.unit_test_data_dir, 'modtran_elevations.txt'),
            ])
        self.setUpTestLinks()

        self.check_files = glob.glob(
            os.path.join(self.unit_test_data_dir, 'LE07*_st_*'))
        self.check_files.extend([os.path.join(self.unit_test_data_dir, 
            'atmospheric_parameters.txt')])
        self.check_files.extend([os.path.join(self.unit_test_data_dir, 
            'used_points.txt')])

    def test_run(self):
        cmd = ["../../src/st_atmospheric_parameters", "--xml", 
                self.xml_filename]
        self.run_test_case(cmd)
        self.compare_updated_xml(".atmos")

class TestBuildSTData(TestST):
    """Test build_st_data.py """

    def setUp(self):
        super(TestBuildSTData, self).setUp()

        self.data_path = os.environ['ST_DATA_DIR']
        self.assertTrue(os.path.exists(self.data_path))

        # This test modifies the XML file, so instead of linking to it, 
        # make a copy.  Also remove it from the link_files list.
        # Use the XML file that contains the ST intermediate outputs.
        shutil.copyfile(
            os.path.join(self.unit_test_data_dir, 
                self.xml_filename + ".build.in"),
            self.xml_filename)
        self.link_files = []

        self.link_files.extend(
            glob.glob(os.path.join(self.unit_test_data_dir, 'LE07*_st_*')))
        self.link_files.extend(
            glob.glob(os.path.join(self.unit_test_data_dir, 'LE07*_emis*img')))
        self.link_files.extend(
            glob.glob(os.path.join(self.unit_test_data_dir, 'LE07*_emis*hdr')))
        self.setUpTestLinks()

        self.check_files.extend(
            glob.glob(os.path.join(self.unit_test_data_dir, 'LE07*_st.*')))

    def test_run(self):
        logging.basicConfig(format=('%(asctime)s.%(msecs)03d %(process)d'
                                    ' %(levelname)-8s'
                                    ' %(filename)s:%(lineno)d:%(funcName)s'
                                    ' -- %(message)s'),
                            datefmt='%Y-%m-%d %H:%M:%S',
                            level=logging.DEBUG)
        current_processor = build_st_data.BuildSTData(
            xml_filename=self.xml_filename)
        current_processor.generate_data()
        self.do_checks()
        self.compare_updated_xml(".st")

if __name__ == '__main__':
    unittest.main(verbosity=2)
