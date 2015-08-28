#------------------------------------------------------------------------------
# Makefile
#
# Simple makefile for building and installing land-surface-temperature.
#------------------------------------------------------------------------------
.PHONY: all install clean

export ESPA_APPLICATION = l5-7_lst

all:
	echo "make all in src..."; \
        (cd src; $(MAKE) -f Makefile)

install: all
	echo "make install in scripts..."; \
        (cd scripts; $(MAKE) -f Makefile install)
	echo "make install in static_data..."; \
        (cd static_data; $(MAKE) -f Makefile install)
	echo "make install in src..."; \
        (cd src; $(MAKE) -f Makefile install)

clean:
	echo "make clean in src..."; \
        (cd clean; $(MAKE) -f Makefile)

