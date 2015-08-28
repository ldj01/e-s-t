#------------------------------------------------------------------------------
# Makefile
#
# Simple makefile for installing the scripts.
#------------------------------------------------------------------------------
.PHONY: all install clean check-environment

export ESPA_APPLICATION = lst_auxillary_data

AUX_INSTALL_PATH = $(PREFIX)/$(ESPA_PROJECT)/$(ESPA_APPLICATION)/bin

all:

install: check-environment
	install -d $(AUX_INSTALL_PATH)
	install -m 755 build_narr_aux_archive_from_CISL_RDA.py $(AUX_INSTALL_PATH)
	install -m 644 example-lst_auxillary.config $(AUX_INSTALL_PATH)

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

