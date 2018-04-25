## Surface Temperature 1.1.0 Release Notes
## Note: The productions produced by this software are considered provisional.
Release Date: May 2018

See git tag [st-rit-v1.1.0]

This project contains application source code for producing Surface Temperature products.

## Product Descriptions
See the [provisional_st_README_V10.pdf](https://edclpdsftp.cr.usgs.gov/downloads/provisional/land_surface_temperature/provisional_st_README_v10.pdf) product guide (which is an unofficial and provisional version) for information about the Surface Temperature products.

## Release Notes
* Version change
* Support global processing (within constraints of available inputs) by using
  MERRA-2 reanalysis input 

## Installation

### Dependencies
* ESPA raw binary libraries, tools, and its dependencies. [Found here](https://github.com/USGS-EROS/espa-product-formatter)
* Python 2.7+ and Numpy/GDAL
* [GDAL](http://www.gdal.org/) 1.11.1
  - The GDAL command line tools are utilized for some of the processing steps.

### Environment Variables
* Required for building this software (For an example see setup-build-environment.sh)
```
export PREFIX="path_to_Installation_Directory"
export XML2INC="path_to_LIBXML2_include_files"
export XML2LIB="path_to_LIBXML2_libraries"
export LZMALIB="path_to_LZMA_libraries_for_linking"
export ZLIBLIB="path_to_ZLIB_libraries_for_linking"
export ESPAINC="path_to_ESPA_PRODUCT_FORMATTER_include_files"
export ESPALIB="path_to_ESPA_PRODUCT_FORMATTER_libraries"
```

### Build Steps
* Clone the repository and replace the defaulted version(master) with this
  version of the software
```
git clone https://github.com/USGS-EROS/espa-surface-temperature.git
cd espa-surface-temperature
git checkout version_<version>
```
* Build and install the application specific software
```
make
make install
```

## Usage
See `surface_temperature.py --help` for command line details.

### Environment Variables
* PATH - May need to be updated to include the following
  - `$PREFIX/bin`
* ST_AUX_DIR - Points to the local NARR data archive.  See [ST Auxiliary Data](st_auxiliary_data/README.md).
  - `export ST_AUX_DIR="/usr/local/auxiliaries/LST/NARR"`
* ST_MERRA_AUX_DIR - Points to the local MERRA-2 data archive.  See [ST Auxiliary Data](st_auxiliary_data/README.md).
  - `export ST_MERRA_AUX_DIR="/usr/local/auxiliaries/LST/MERRA2"`
* ST_DATA_DIR - Points to the installed static file
  - `export ST_DATA_DIR="/usr/local/espa-surface-temperature/st/static_data"`
* MODTRAN_PATH - Points to the installed MODTRAN location
  - `export MODTRAN_PATH="/usr/local/bin"`
* MODTRAN_DATA_DIR - Points to the directory containing the MODTRAN "DATA" directory
  - `export MODTRAN_DATA_DIR="/usr/local/auxiliaries/MODTRAN_DATA"`
* ASTER_GED_SERVER_NAME
  - `export ASTER_GED_SERVER_NAME="e4ftl01.cr.usgs.gov"`
* ASTER_GED_SERVER_PATH
  - `export ASTER_GED_SERVER_PATH="/ASTT/AG100.003/2000.01.01/"`

### Data Processing Requirements
This version of the Surface Temperature application requires the input products to be in the ESPA internal file format.

The following input data are required to generate the Surface Temperature product:
* Level 1 Product
  - The thermal band from the level 1 product is used.
* Top of Atmosphere Reflectance (TOA)
  - TOA products can be generated using the [LEDAPS](https://github.com/USGS-EROS/espa-surface-reflectance) or [L8_SR](https://github.com/USGS-EROS/espa-surface-reflectance) software found in our [espa-surface-reflectance](https://github.com/USGS-EROS/espa-surface-reflectance) project.  TOA products can also be generated through our on-demand processing system [ESPA](https://espa.cr.usgs.gov).  Be sure to select the ENVI output format if using our [ESPA](https://espa.cr.usgs.gov) on-demand system.
* Elevation
  - Elevation data is expected to be in the same projection, resolution, and image dimensions as the TOA products as well as identified in the XML metadata file and in the ENVI file format.  The filename is also expected be `<sceneid>_elevation.img`: where <b>sceneid</b> is the Landsat Scene ID of the data being processed.
* ASTER GED
  - ASTER GED data can be [found here](https://lpdaac.usgs.gov/data_access/data_pool).  External users will need to set up an alternate ASTER GED data retrieval method.  One potential method is to create a NASA Earthdata account as [described here] (http://e4ftl01.cr.usgs.gov/ASTT), and update the automated ST procedure to supply the login information before accessing the ASTER GED data.
* North American Regional Reanalysis (NARR)
  - For NARR data, it would be best to use the `st_auxiliary_data` software provided in this project to download and build your own archive for the dates you require.  This software archives a reduced set of parameters from each source file, only using the parameters required for ST generation.  See [ST Auxiliary Data](st_auxiliary_data/README.md).
* Modern-Era Retrospective analysis for Research and Applications, Version 2 
  (MERRA-2)
  - For MERRA-2 data, it would be best to use the `st_auxiliary_data` software provided in this project to download and build your own archive for the dates you require.  This software archives a reduced set of parameters from each source file, only using the parameters required for ST generation.  See [ST Auxiliary Data](st_auxiliary_data/README.md).

### Data Postprocessing
After compiling the [espa-product-formatter](https://github.com/USGS-EROS/espa-product-formatter) libraries and tools, the `convert_espa_to_gtif` and `convert_espa_to_hdf` command-line tools can be used to convert the ESPA internal file format to HDF or GeoTIFF.  Otherwise the data will remain in the ESPA internal file format, which includes each band in the ENVI file format (i.e. raw binary file with associated ENVI header file) and an overall XML metadata file.

### Support Information
This project is unsupported software provided by the U.S. Geological Survey (USGS) Earth Resources Observation and Science (EROS) Land Satellite Data Systems (LSDS) Project. For questions regarding products produced by this source code, please contact us at custserv@usgs.gov. 

### Disclaimer
This software is preliminary or provisional and is subject to revision. It is being provided to meet the need for timely best science. The software has not received final approval by the U.S. Geological Survey (USGS). No warranty, expressed or implied, is made by the USGS or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. The software is provided on the condition that neither the USGS nor the U.S. Government shall be held liable for any damages resulting from the authorized or unauthorized use of the software.
