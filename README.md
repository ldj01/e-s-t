
This project implements algorithms for determing Land Surface Temperature.

Algorithms for specific sensors are located in sensor specific sub-directories.

See the sensor specific sub-directories for more details and usage examples.


Installation Instructions:
    - All sensor applications can be built and installed from the top-level
      directory, or at the users choice only a single individual sensors
      applications could be built and installed.

    - To build all sensor applications.
        make all

    - To install all sensor applications.
        make install

    - To build a single sensor application
        make <sensor>

    - To install a single sensor application
        make <sensor>-install

    Available sensor options are:

        l5-7_lst


This project is hosted by the US Geological Survey (USGS) Earth Resources Observation and Science (EROS) Land Satellite Data Systems (LSDS) Science Research and Development (LSRD) Project. For questions regarding this source code, please contact the Landsat Contact Us page and specify USGS CDR/ECV in the "Regarding" section. https://landsat.usgs.gov/contactus.php
