#------------------------------------------------------------------------------
# Makefile
#
# Simple makefile for installing the scripts.
#------------------------------------------------------------------------------
.PHONY: all install clean

# Inherit from upper-level make.config
TOP = ../..
include $(TOP)/make.config

script_install_path = $(l57_algorithm_dir)/bin

SCRIPT_EXECUTABLES = \
    lst_core_processing.py \
    extract_modtran_results.py

SCRIPT_IMPORTS = \
    build_lst_data.py \
    estimate_landsat_emissivity.py \
    lst_utilities.py

all:

install:
	install -d $(link_path)
	install -d $(script_install_path)
	@for file in $(SCRIPT_EXECUTABLES); do \
            echo "  installing $$file"; \
            install -m 755 $$file $(script_install_path); \
            ln -sf ../$(project_name)/$(l57_algorithm)/bin/$$file $(link_path)/$$file; \
        done;
	@for file in $(SCRIPT_IMPORTS); do \
            echo "  installing $$file"; \
            install -m 644 $$file $(script_install_path); \
            ln -sf ../$(project_name)/$(l57_algorithm)/bin/$$file $(link_path)/$$file; \
        done;

clean:

