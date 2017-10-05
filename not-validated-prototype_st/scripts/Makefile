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
    st_determine_grid_points.py \
    st_build_modtran_input.py \
    st_extract_auxiliary_narr_data.py \
    st_generate_distance_to_cloud.py \
    st_generate_products.py \
    st_generate_qa.py \
    st_run_modtran.py \
    estimate_landsat_emissivity.py

SCRIPT_IMPORTS = \
    build_st_data.py \
    st_exceptions.py \
    st_grid_points.py \
    st_utilities.py

#-----------------------------------------------------------------------------
all:

#-----------------------------------------------------------------------------
install:
	install -d $(link_path)
	install -d $(st_install_path)
	@for file in $(SCRIPTS); do \
            echo "  installing $$file"; \
            install -m 755 $$file $(st_install_path) || exit 1; \
            ln -sf $(st_link_source_path)/$$file $(link_path)/$$file; \
        done;
	@for file in $(SCRIPT_IMPORTS); do \
            echo "  installing $$file"; \
            install -m 644 $$file $(st_install_path) || exit 1; \
            ln -sf $(st_link_source_path)/$$file $(link_path)/$$file; \
        done;

#-----------------------------------------------------------------------------
clean:

#-----------------------------------------------------------------------------
check:
	@cd unit-tests; python unit-tests.py

