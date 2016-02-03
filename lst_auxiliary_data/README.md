## Land Surface Temperature Auxiliary Archive Generation

The scripts in this directory facilitate building an archive of NARR data used by the core Land Surface Temperature processing during product generation.

#### build_narr_aux_archive_from_CISL_RDA.py

This script is used to retrieve data from the NARR archive located at http://rda.ucar.edu   This data is provided in 1 to 4 day increments depending on the year and days in the month.  Most files are 3day files.  When a file is downloaded and processed, all days with the file will be processed and added(or will update) the local archive.  See http://rda.ucar.edu for more details about the contents of the data.

#### update_narr_aux_data.py

This script is used to update the archive on a daily basis from http://ftp.cpc.ncep.noaa.gov/wd51we/NARR_archive where files are available for the current year in 3 hour increments.  However there is a delay before new files are provided as time is needed to wait for input data and to generate the parameters.  See http://rda.ucar.edu for more details about the contents of the data.
