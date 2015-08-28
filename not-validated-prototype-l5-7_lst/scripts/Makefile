#------------------------------------------------------------------------------
# Makefile
#
# Simple makefile for installing the scripts.
#------------------------------------------------------------------------------
.PHONY: all install clean check-environment

SCRIPT_INSTALL_PATH = $(PREFIX)/$(ESPA_PROJECT)/$(ESPA_APPLICATION)/bin

SCRIPT_FILES = \
    l5_7_lst.py \
    l5_7_build_lst_data.py \
    l5_7_estimate_landsat_emissivity.py \
    lst_extract_modtran_results.py \
    lst_utilities.py \

all:

install: check-environment
	install -d $(SCRIPT_INSTALL_PATH)
	@for file in $(SCRIPT_FILES); do \
            echo "  installing $$file"; \
            install -m 755 $$file $(SCRIPT_INSTALL_PATH); \
        done;

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

