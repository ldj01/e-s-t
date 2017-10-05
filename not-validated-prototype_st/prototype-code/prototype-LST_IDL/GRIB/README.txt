To grab the NARR data (in the .grib format) for a particular date/time/variable, use a text editor to edit the file "script" and then type ./script.
The data will be downloaded to a file named out.grb

To convert from .grib to .txt, type the following...
./a.out out.grb -d 1 -text -o output_filename.txt

