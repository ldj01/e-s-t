## NOTE - The products produced by this software have not been validated and are considered prototype.

## Land Surface Temperature 0.1.0 Release Notes
Release Date: January 2017

See git tag [2017_Jan]

This project contains application source code for producing Land Surface Temperature products.

## Product Descriptions
See the [lst_readme_vX.X_for_evaluators.pdf](http://espa.cr.usgs.gov/downloads/provisional/land_surface_temperature) product guide (which is an unofficial and provisional version) for information about the Land Surface Temperature products.

## Release Notes
* Version change
* The software was significantly restructured to only process the needed 
  points through MODTRAN 
* Some modules were converted from C to Python  
* Some modules were consolidated into fewer files
* C module memory handling was improved using suggestions from Valgrind 
* Python module formatting was updated using suggestions from Pylint

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
git clone https://github.com/USGS-EROS/espa-land-surface-temperature.git
cd espa-land-surface-temperature
git checkout version_<version>
```
* Build and install the application specific software
```
make
make install
```

## Usage
See `land_surface_temperature.py --help` for command line details.

### Environment Variables
* PATH - May need to be updated to include the following
  - `$PREFIX/bin`
* LST_AUX_DIR - Points to our NARR data archive
  - `export LST_AUX_DIR="/usr/local/auxiliaries/LST/NARR"`
* LST_DATA_DIR - Points to the installed static file
  - `export LST_DATA_DIR="/usr/local/espa-land-surface-temperature/lst/static_data"`
* MODTRAN_PATH - Points to the installed MODTRAN location
  - `export MODTRAN_PATH="/usr/local/bin"`
* MODTRAN_DATA_DIR - Points to the directory containing the MODTRAN "DATA" directory
  - `export MODTRAN_DATA_DIR="/usr/local/auxiliaries/MODTRAN_DATA"`
* ASTER_GED_SERVER_NAME
  - `export ASTER_GED_SERVER_NAME="e4ftl01.cr.usgs.gov"`
* ASTER_GED_SERVER_PATH
  - `export ASTER_GED_SERVER_PATH="/ASTT/AG100.003/2000.01.01/"`

### Data Processing Requirements
This version of the Land Surface Temperature application requires the input products to be in the ESPA internal file format.

The following input data are required to generate the Land Surface Temperature product:
* Top of Atmosphere Reflectance (TOA)
  - TOA products can be generated using the [LEDAPS](https://github.com/USGS-EROS/espa-surface-reflectance) or [L8_SR](https://github.com/USGS-EROS/espa-surface-reflectance) software found in our [espa-surface-reflectance](https://github.com/USGS-EROS/espa-surface-reflectance) project.  Or through our ondemand processing system [ESPA](https://espa.cr.usgs.gov), be sure to select the ENVI output format.
* Elevation
  - Elevation data is expected to be in the same projection and image dimensions as the TOA products and in the ENVI file format.  The filename is also expected be `<sceneid>_dem.img`: where <b>sceneid</b> is the Landsat Scene ID of the data being processed.
* ASTER GED
  - ASTER GED data can be [found here](https://lpdaac.usgs.gov/data_access/data_pool).  However it will automatically be retrieved for you as needed and cleaned up.
* North American Regional Reanalysis (NARR)
  - For NARR data, it would be best to utilize the `lst_aux_data` software provided in this project to download and build your own archive for the dates you require.

### Data Postprocessing
After compiling the [espa-product-formatter](https://github.com/USGS-EROS/espa-product-formatter) libraries and tools, the `convert_espa_to_gtif` and `convert_espa_to_hdf` command-line tools can be used to convert the ESPA internal file format to HDF or GeoTIFF.  Otherwise the data will remain in the ESPA internal file format, which includes each band in the ENVI file format (i.e. raw binary file with associated ENVI header file) and an overall XML metadata file.

## More Information
This project is provided by the US Geological Survey (USGS) Earth Resources Observation and Science (EROS) Land Satellite Data Systems (LSDS) Science Research and Development (LSRD) Project. For questions regarding products produced by this source code, please contact the Landsat Contact Us page and specify USGS CDR/ECV in the "Regarding" section. https://landsat.usgs.gov/contactus.php
