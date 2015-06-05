#------------------------------------------------------------------------------
# Makefile
#
# Simple makefile for building and installing land-surface-temperature
# applications.
#------------------------------------------------------------------------------

all: l5-7_lst
	echo "make all in scripts"; \
        (cd scripts; $(MAKE) all -f Makefile);

install: l5-7_lst-install
	echo "make install in scripts"; \
        (cd scripts; $(MAKE) install -f Makefile);

clean: l5-7_lst-clean
	echo "make clean in scripts"; \
        (cd scripts; $(MAKE) clean -f Makefile);

l5-7_lst: l5-7_lst-all

l5-7_lst-all:
	echo "make all in not-validated-prototype-l5-7_lst"; \
        (cd not-validated-prototype-l5-7_lst; \
        $(MAKE) all -f Makefile);

l5-7_lst-install:
	echo "make install in not-validated-prototype-l5-7_lst"; \
        (cd not-validated-prototype-l5-7_lst; \
        $(MAKE) install -f Makefile);

l5-7_lst-clean:
	echo "make clean in not-validated-prototype-l5-7_lst"; \
        (cd not-validated-prototype-l5-7_lst; \
        $(MAKE) clean -f Makefile);
