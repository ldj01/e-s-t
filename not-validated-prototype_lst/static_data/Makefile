#-----------------------------------------------------------------------------
# Makefile
#
# Simple makefile for installing the static data.
#-----------------------------------------------------------------------------
.PHONY: all install clean

# Inherit from upper-level make.config
TOP = ../..
include $(TOP)/make.config

static_install_path = $(lst_algorithm_dir)/static_data

STATIC_DATA_FILES = \
    modtran_head.txt                  \
    modtran_tail.txt                  \
    narr_coordinates.txt              \
    std_mid_lat_summer_atmos.txt      \
    L4_Brightness_Temperature_LUT.txt \
    L4_Spectral_Response.txt          \
    L5_Brightness_Temperature_LUT.txt \
    L5_Spectral_Response.txt          \
    L7_Brightness_Temperature_LUT.txt \
    L7_Spectral_Response.txt          \
    L8_Brightness_Temperature_LUT.txt \
    L8_Spectral_Response.txt

all:

install:
	install -d $(static_install_path)
	@for file in $(STATIC_DATA_FILES); do \
            echo "  installing $$file"; \
            install -m 644 $$file $(static_install_path) || exit 1; \
        done;

clean:

