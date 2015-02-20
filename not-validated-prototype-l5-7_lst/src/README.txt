CONFIGURATION:
1. Need $LST_DATA environmental variable be set in which all needed input text
   files are located.
2. Need $MODTRAN_PATH be set in which MODTRAN executable file is located and
   DATA subdiectory exists.
3. Need $BIN directory be set with scene_based_lst executable file and all
   needed scripts be copied over.
4. wgrib needs to be in the path.


PROCESSING:
1. Unavailable strings
    "(Satellite not supported with LST processing)"
    "(Sensor not supported with LST processing)"


QUESTIONS:
1. Should we be using Radiance or TOA for the thermal data?
    a. If Radiance, why is 0.044 added for L5 thermal only??????

