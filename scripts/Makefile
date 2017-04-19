#-----------------------------------------------------------------------------
# Makefile
#
# Simple makefile for installing the scripts.
#-----------------------------------------------------------------------------
.PHONY: all install clean

# Inherit from upper-level make.config
TOP = ..
include $(TOP)/make.config

# Define the script installation path
script_install_path = $(espa_project_dir)/bin

all:

install:
	install -d $(link_path)
	install -d $(script_install_path)
	install -m 755 land_surface_temperature.py $(script_install_path) || exit 1
	ln -sf ../$(project_name)/bin/land_surface_temperature.py $(link_path)/land_surface_temperature.py

clean:

