## NOTE - Since we are still in the very early stages of development, the information contained in this README may not be correct, however it will be updated from time to time.

## NOTE - The products produced by this software have not been validated and are considered prototype.

## Land Surface Temperature 0.0.1 Release Notes
Release Date: MONTH DAY, YEAR

This application produces Land Surface Temperature products for Landsats 5, and 7.

## Product Descriptions
See the [PLACE HOLDER](butter) product guide for information about the Land Surface Temeprature products.

## Release Notes
* Initial implementation

## Installation

### Dependencies
* ESPA raw binary libraries, tools, and it's dependencies. [Found here](https://github.com/USGS-EROS/espa-product-formatter)
* Python 2.7+ and Scipy/Numpy/GDAL
* [GDAL](http://www.gdal.org/) 1.11.1
  - The command line tools are utilized for some of the processing steps.

### Environment Variables
* Required for building this software
```
export PREFIX="path_to_Installation_Directory"
export XML2INC="path_to_LIBXML2_include_files"
export XML2LIB="path_to_LIBXML2_libraries"
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
* ASTER_GED_SERVER_NAME
  - `export ASTER_GED_SERVER_NAME="e4ftl01.cr.usgs.gov"`
* ASTER_GED_SERVER_PATH
  - `export ASTER_GED_SERVER_PATH="/ASTT/AG100.003/2000.01.01/"`

### Data Processing Requirements
This version of the Land Surface Temperature application requires the input products to be in the ESPA internal file format.

The following input data are required to generate the spectral indicies products:
* Top of Atmosphere Reflectance (TOA)
  - TOA products can be generated using the [LEDAPS](https://github.com/USGS-EROS/espa-surface-reflectance) or [L8_SR](https://github.com/USGS-EROS/espa-surface-reflectance) software found in our [espa-surface-reflectance](https://github.com/USGS-EROS/espa-surface-reflectance) project.  Or through our ondemand processing system [ESPA](https://espa.cr.usgs.gov), be sure to select the ENVI output format.
* Elevation
  - Elevation data is expected to be in the same projection and image dimensions as the TOA products and in the ENVI file format.  The filename is also expected be `<sceneid>_dem.img`: where <b>sceneid</b> is the Landsat Scene ID of the data being processed.
* ASTER GED
  - ASTER GED data can be [found here](https://lpdaac.usgs.gov/data_access/data_pool).  However it will automatically be retrieved for you as needed and cleaned up.
* North American Regional Reanalysis (NARR)
  - NARR data, it would be best to utilize the `lst_aux_data` software provided in this project to download and build your own archive for the dates you require.

### Data Postprocessing
After compiling the [espa-product-formatter](https://github.com/USGS-EROS/espa-product-formatter) libraries and tools, the `convert_espa_to_gtif` and `convert_espa_to_hdf` command-line tools can be used to convert the ESPA internal file format to HDF or GeoTIFF.  Otherwise the data will remain in the ESPA internal file format, which includes each band in the ENVI file format (i.e. raw binary file with associated ENVI header file) and an overall XML metadata file.

## More Information
This project is provided by the US Geological Survey (USGS) Earth Resources Observation and Science (EROS) Land Satellite Data Systems (LSDS) Science Research and Development (LSRD) Project. For questions regarding products produced by this source code, please contact the Landsat Contact Us page and specify USGS CDR/ECV in the "Regarding" section. https://landsat.usgs.gov/contactus.php
