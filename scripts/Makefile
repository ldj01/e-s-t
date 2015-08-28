#------------------------------------------------------------------------------
# Makefile
#
# Simple makefile for installing the scripts.
#------------------------------------------------------------------------------
.PHONY: all install clean check-environment

SCRIPT_INSTALL_PATH = $(PREFIX)/$(ESPA_PROJECT)/bin

all:

install: check-environment
	install -d $(SCRIPT_INSTALL_PATH)
	install -m 755 land_surface_temperature.py $(SCRIPT_INSTALL_PATH)

clean:

check-environment:
ifndef PREFIX
    $(error PREFIX is undefined)
endif
ifndef ESPA_PROJECT
    $(error ESPA_PROJECT is undefined)
endif

