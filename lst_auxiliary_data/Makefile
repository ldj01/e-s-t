#-----------------------------------------------------------------------------
# Makefile
#
# Simple makefile for installing the scripts.
#-----------------------------------------------------------------------------
.PHONY: all install clean

# Inherit from upper-level make.config
TOP = ..
include $(TOP)/make.config

SCRIPTS = lst_aux_narr_from_CISL_RDA_archive.py \
          lst_aux_update_narr_data.py
MODULES = lst_aux_utilities.py \
          lst_aux_config.py \
          lst_aux_exception.py \
          lst_aux_parameters.py \
          lst_aux_http_session.py \
          lst_aux_version.py \
          lst_aux_logging.py


#-----------------------------------------------------------------------------
all:

#-----------------------------------------------------------------------------
install:
	install -d $(link_path)
	install -d $(lst_auxiliary_install_path)
	@for script in $(SCRIPTS); do \
            cmd="install -m 755 $$script $(lst_auxiliary_install_path)"; \
            echo "$$cmd"; $$cmd; \
            cmd="ln -sf $(lst_auxiliary_link_source_path)/$$script $(link_path)/$$script"; \
            echo "$$cmd"; $$cmd; \
        done
	@for module in $(MODULES); do \
            cmd="install -m 644 $$module $(lst_auxiliary_install_path)"; \
            echo "$$cmd"; $$cmd; \
            cmd="ln -sf $(lst_auxiliary_link_source_path)/$$module $(link_path)/$$module"; \
            echo "$$cmd"; $$cmd; \
        done

#-----------------------------------------------------------------------------
clean:

