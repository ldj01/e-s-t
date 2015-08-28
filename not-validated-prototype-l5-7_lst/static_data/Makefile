#------------------------------------------------------------------------------
# Makefile
#
# Simple makefile for installing the static data.
#------------------------------------------------------------------------------
.PHONY: check-environment

STATIC_INSTALL_PATH = $(PREFIX)/$(ESPA_PROJECT)/$(ESPA_APPLICATION)/static_data

STATIC_DATA_FILES = \
    modtran_head.txt                  \
    modtran_tail.txt                  \
    narr_coordinates.txt              \
    std_mid_lat_summer_atmos.txt      \
    L5_Brightness_Temperature_LUT.txt \
    L5_Spectral_Response.txt          \
    L7_Brightness_Temperature_LUT.txt \
    L7_Spectral_Response.txt          \
    L8_Spectral_Response.txt

all:

install: check-environment
	install -d $(STATIC_INSTALL_PATH)
	@for file in $(STATIC_DATA_FILES); do \
            echo "  installing $$file"; \
            install -m 644 $$file $(STATIC_INSTALL_PATH); \
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

