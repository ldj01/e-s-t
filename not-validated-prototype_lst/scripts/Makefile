#-----------------------------------------------------------------------------
# Makefile
#
# Simple makefile for installing the scripts.
#-----------------------------------------------------------------------------
.PHONY: all install clean check

# Inherit from upper-level make.config
TOP = ../..
include $(TOP)/make.config

SCRIPTS = \
    lst_core_processing.py \
    lst_determine_grid_points.py \
    extract_modtran_results.py

SCRIPT_IMPORTS = \
    build_lst_data.py \
    estimate_landsat_emissivity.py \
    extract_auxiliary_narr_data.py \
    lst_environment.py \
    lst_exceptions.py \
    lst_utilities.py

#-----------------------------------------------------------------------------
all:

#-----------------------------------------------------------------------------
install:
	install -d $(link_path)
	install -d $(lst_install_path)
	@for file in $(SCRIPTS); do \
            echo "  installing $$file"; \
            install -m 755 $$file $(lst_install_path); \
            ln -sf $(lst_link_source_path)/$$file $(link_path)/$$file; \
        done;
	@for file in $(SCRIPT_IMPORTS); do \
            echo "  installing $$file"; \
            install -m 644 $$file $(lst_install_path); \
            ln -sf $(lst_link_source_path)/$$file $(link_path)/$$file; \
        done;

#-----------------------------------------------------------------------------
clean:

#-----------------------------------------------------------------------------
check:
	@cd unit-tests; python unit-tests.py

