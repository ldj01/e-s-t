#-----------------------------------------------------------------------------
# Makefile
#
# Simple makefile for building and installing land-surface-temperature
# applications.
#-----------------------------------------------------------------------------
.PHONY: check-environment all install clean all-script install-script clean-script all-l5-7 install-l5-7 clean-l5-7 install-aux

include make.config

DIR_L5-7 = not-validated-prototype-l5-7_lst
DIR_AUX = lst_auxiliary_data

all: all-script all-l5-7

install: check-environment install-script install-l5-7

clean: clean-script clean-l5-7

#------------------------------------------------------------------------------
all-script:
	echo "make all in scripts"; \
        (cd scripts; $(MAKE) all);

install-script: check-environment
	echo "make install in scripts"; \
        (cd scripts; $(MAKE) install);

clean-script:
	echo "make clean in scripts"; \
        (cd scripts; $(MAKE) clean);

#------------------------------------------------------------------------------
all-l5-7: all-script
	echo "make all in not-validated-prototype-l5-7_lst"; \
        (cd $(DIR_L5-7); $(MAKE) all);

install-l5-7: check-environment install-script
	echo "make install in not-validated-prototype-l5-7_lst"; \
        (cd $(DIR_L5-7); $(MAKE) install);

clean-l5-7: clean-script
	echo "make clean in not-validated-prototype-l5-7_lst"; \
        (cd $(DIR_L5-7); $(MAKE) clean);

#------------------------------------------------------------------------------
install-aux:
	echo "make install in lst_auxiliary_data"; \
        (cd $(DIR_AUX); $(MAKE) install);

check-environment:
ifndef PREFIX
    $(error Environment variable PREFIX is not defined)
endif

