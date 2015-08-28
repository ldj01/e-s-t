## Land Surface Temperature Auxiliary Archive Generation

The scripts in this directory facilitate building an archive of NARR data used by the core Land Surface Temperature processing during product generation.

#### build_narr_aux_archive_from_CISL_RDA.py

This script is used to retrieve data from the NARR archive located at http://rda.ucar.edu   This data is provided in 1 to 4 day increments depending on the year and days in the month.  Most files are 3day files.  When a file is downloaded and processed, all days with the file will be processed and added(or will update) the local archive.

#### ??????????.py
