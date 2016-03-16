#-----------------------------------------------------------------------------
# Makefile
#
# Simple makefile for building and installing land-surface-temperature
# applications.
#-----------------------------------------------------------------------------
.PHONY: check-environment all install clean all-script install-script clean-script all-lst install-lst clean-lst all-aux install-aux

include make.config

DIR_LST = not-validated-prototype_lst
DIR_AUX = lst_auxiliary_data

all: all-script all-lst

install: check-environment install-script install-lst

clean: clean-script clean-lst clean-aux

#-----------------------------------------------------------------------------
all-script:
	echo "make all in scripts"; \
        (cd scripts; $(MAKE) all);

install-script: check-environment
	echo "make install in scripts"; \
        (cd scripts; $(MAKE) install);

clean-script:
	echo "make clean in scripts"; \
        (cd scripts; $(MAKE) clean);

#-----------------------------------------------------------------------------
all-lst: all-script
	echo "make all in not-validated-prototype_lst"; \
        (cd $(DIR_LST); $(MAKE) all);

install-lst: check-environment install-script
	echo "make install in not-validated-prototype_lst"; \
        (cd $(DIR_LST); $(MAKE) install);

clean-lst: clean-script
	echo "make clean in not-validated-prototype_lst"; \
        (cd $(DIR_LST); $(MAKE) clean);

#-----------------------------------------------------------------------------
all-aux:
	echo "make install in lst_auxiliary_data"; \
        (cd $(DIR_AUX); $(MAKE));

install-aux:
	echo "make install in lst_auxiliary_data"; \
        (cd $(DIR_AUX); $(MAKE) install);

clean-aux:
	echo "make install in lst_auxiliary_data"; \
        (cd $(DIR_AUX); $(MAKE) clean);

#-----------------------------------------------------------------------------
check-environment:
ifndef PREFIX
    $(error Environment variable PREFIX is not defined)
endif

