#------------------------------------------------------------------------------
# Makefile
#
# Simple makefile for installing the scripts.
#------------------------------------------------------------------------------

# Target for the executable
all:

install:
	install -d $(PREFIX)/bin
	install -m 755 l5-7_lst.py $(PREFIX)/bin
	install -m 755 l5-7_build_lst_data.py $(PREFIX)/bin
	install -m 755 download_narr_aux_data.py $(PREFIX)/bin
	install -m 755 lst_extract_modtran_results.py $(PREFIX)/bin
	install -m 755 lst_utilities.py $(PREFIX)/bin
	install -m 755 l5-7_landsat_emissivity_from_aster_ged.py $(PREFIX)/bin

clean:

