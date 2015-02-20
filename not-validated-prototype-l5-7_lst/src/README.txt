CONFIGURATION:
1. New environment variables:

       $LST_DATA_DIR - Where all the static files are located.

       $LST_AUX_DIR  - Where the input grib files are archived.

       $MODTRAN_PATH - Where the MODTRAN executable file is located and where
                       the "DATA" subdiectory exists.


2. Executables needed in the $PATH:
   lst
   do_lst.py
   wgrib


3. New tools to be added to the system.
   MODTRAN x.x.x
   wgrib


PROCESSING:
1. Unavailable strings
    "(Satellite not supported with LST processing)"
    "(Sensor not supported with LST processing)"


QUESTIONS:
1. Should we be using Radiance or TOA for the thermal data?
    a. If Radiance, why is 0.044 added for L5 thermal only??????
    b. If TOA, probably don't need the 0.044 to be added??????

