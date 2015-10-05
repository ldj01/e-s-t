#------------------------------------------------------------------------------
# Makefile
#
# Simple makefile for installing the scripts.
#------------------------------------------------------------------------------
.PHONY: all install clean

# Inherit from upper-level make.config
TOP = ..
include $(TOP)/make.config

aux_install_path = $(lst_auxiliary_dir)/bin

all:

install:
	install -d $(link_path)
	install -d $(aux_install_path)
	install -m 755 build_narr_aux_archive_from_CISL_RDA.py $(aux_install_path)
	install -m 755 update_narr_aux_data.py $(aux_install_path)
	install -m 644 example-lst_auxiliary.config $(aux_install_path)
	ln -sf ../$(project_name)/$(lst_auxiliary)bin/build_narr_aux_archive_from_CISL_RDA.py $(link_path)/build_narr_aux_archive_from_CISL_RDA.py
	ln -sf ../$(project_name)/$(lst_auxiliary)bin/update_narr_aux_data.py $(link_path)/update_narr_aux_data.py
	ln -sf ../$(project_name)/$(lst_auxiliary)bin/example-lst_auxiliary.config $(link_path)/example-lst_auxiliary.config

clean:

