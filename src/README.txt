CONFIGURATION:
1. New environment variables:

       $LST_DATA_DIR - Where all the static files are located.

       $LST_AUX_DIR  - Where the input grib files are archived.

       $MODTRAN_PATH - Where the MODTRAN executable file is located.

       $MODTRAN_DATA_DIR - Exactly the "path/to/DATA" diectory.


2. Executables needed in the $PATH:
   lst
   do_lst.py
   lst_extract_modtran_results.py


3. New tools to be added to the system and path.
   MODTRAN x.x.x
   wgrib


PROCESSING:
1. Unavailable strings
    "(Satellite not supported with LST processing)"
    "(Sensor not supported with LST processing)"


QUESTIONS:
1. Why is 0.044 added for L5 thermal only??????

