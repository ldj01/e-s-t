#------------------------------------------------------------------------------
# Makefile
#
# Simple makefile for installing the static data.
#------------------------------------------------------------------------------
.PHONY: check-environment

STATIC_INSTALL_PATH = $(PREFIX)/$(ESPA_PROJECT)/$(ESPA_APPLICATION)/static_data

all:

install: check-environment
	install -d $(STATIC_INSTALL_PATH)
	install -m 644 modtran_head.txt $(STATIC_INSTALL_PATH)
	install -m 644 modtran_tail.txt $(STATIC_INSTALL_PATH)
	install -m 644 narr_coordinates.txt $(STATIC_INSTALL_PATH)
	install -m 644 std_mid_lat_summer_atmos.txt $(STATIC_INSTALL_PATH)
	install -m 644 L5_Brightness_Temperature_LUT.txt $(STATIC_INSTALL_PATH)
	install -m 644 L5_Spectral_Response.txt $(STATIC_INSTALL_PATH)
	install -m 644 L7_Brightness_Temperature_LUT.txt $(STATIC_INSTALL_PATH)
	install -m 644 L7_Spectral_Response.txt $(STATIC_INSTALL_PATH)
	install -m 644 L8_Spectral_Response.txt $(STATIC_INSTALL_PATH)

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

