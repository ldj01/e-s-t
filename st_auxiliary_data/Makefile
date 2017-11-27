#-----------------------------------------------------------------------------
# Makefile
#
# Simple makefile for installing the scripts.
#-----------------------------------------------------------------------------
.PHONY: all install clean

# Inherit from upper-level make.config
TOP = ..
include $(TOP)/make.config

SCRIPTS = st_aux_narr_from_CISL_RDA_archive.py \
          st_aux_update_narr_data.py
MODULES = st_aux_utilities.py \
          st_aux_config.py \
          st_aux_exception.py \
          st_aux_parameters.py \
          st_aux_http_session.py \
          st_aux_version.py \
          st_aux_logging.py


#-----------------------------------------------------------------------------
all:

#-----------------------------------------------------------------------------
install:
	install -d $(link_path)
	install -d $(st_auxiliary_install_path)
	@for script in $(SCRIPTS); do \
            cmd="install -m 755 $$script $(st_auxiliary_install_path)"; \
            echo "$$cmd"; $$cmd || exit 1; \
            cmd="ln -sf $(st_auxiliary_link_source_path)/$$script $(link_path)/$$script"; \
            echo "$$cmd"; $$cmd; \
        done
	@for module in $(MODULES); do \
            cmd="install -m 644 $$module $(st_auxiliary_install_path)"; \
            echo "$$cmd"; $$cmd || exit 1; \
            cmd="ln -sf $(st_auxiliary_link_source_path)/$$module $(link_path)/$$module"; \
            echo "$$cmd"; $$cmd; \
        done

#-----------------------------------------------------------------------------
clean:

