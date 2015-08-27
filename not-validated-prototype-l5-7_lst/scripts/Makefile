#------------------------------------------------------------------------------
# Makefile
#
# Simple makefile for installing the scripts.
#------------------------------------------------------------------------------
.PHONY: check-environment

SCRIPT_INSTALL_PATH = $(PREFIX)/$(ESPA_PROJECT)/$(ESPA_APPLICATION)/bin

# Target for the executable
all:

install: check-environment
	install -d $(SCRIPT_INSTALL_PATH)
	install -m 755 l5_7_lst.py $(SCRIPT_INSTALL_PATH)
	install -m 755 l5_7_build_lst_data.py $(SCRIPT_INSTALL_PATH)
	install -m 755 l5_7_estimate_landsat_emissivity.py $(SCRIPT_INSTALL_PATH)
	install -m 755 lst_extract_modtran_results.py $(SCRIPT_INSTALL_PATH)
	install -m 644 lst_utilities.py $(SCRIPT_INSTALL_PATH)

clean:

check-environment:
ifndef PREFIX
    $(error PREFIX is undefined)
endif
ifndef ESPA_PROJECT
    $(error ESPA_PROJECT is undefined)
endif
ifndef ESPA_APPLICATION
    $(error ESPA_APPLICATION is undefined)
endif

