#------------------------------------------------------------------------------
# Makefile
#
# Simple makefile for building and installing land-surface-temperature
# applications.
#------------------------------------------------------------------------------
.PHONY: check-environment all install clean all-script install-script clean-script all-l5-7 install-l5-7 clean-l5-7 install-aux

include make.config

MAKEFILE_NAME = Makefile

DIR_L5-7 = not-validated-prototype-l5-7_lst
DIR_AUX = lst_auxiliary_data

all: all-script all-l5-7

install: check-environment install-script install-l5-7

clean: clean-script clean-l5-7

#------------------------------------------------------------------------------
all-script:
	echo "make all in scripts"; \
        (cd scripts; $(MAKE) all -f $(MAKEFILE_NAME));

install-script: check-environment
	echo "make install in scripts"; \
        (cd scripts; $(MAKE) install -f $(MAKEFILE_NAME));

clean-script:
	echo "make clean in scripts"; \
        (cd scripts; $(MAKE) clean -f $(MAKEFILE_NAME));

#------------------------------------------------------------------------------
all-l5-7: all-script
	echo "make all in not-validated-prototype-l5-7_lst"; \
        (cd $(DIR_L5-7); $(MAKE) all -f $(MAKEFILE_NAME));

install-l5-7: check-environment install-script
	echo "make install in not-validated-prototype-l5-7_lst"; \
        (cd $(DIR_L5-7); $(MAKE) install -f $(MAKEFILE_NAME));

clean-l5-7: clean-script
	echo "make clean in not-validated-prototype-l5-7_lst"; \
        (cd $(DIR_L5-7); $(MAKE) clean -f $(MAKEFILE_NAME));

#------------------------------------------------------------------------------
install-aux:
	echo "make install in lst_auxiliary_data"; \
        (cd $(DIR_AUX); $(MAKE) install -f $(MAKEFILE_NAME));

check-environment:
ifndef PREFIX
    $(error Environment variable PREFIX is not defined)
endif

