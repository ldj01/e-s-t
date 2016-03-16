#-----------------------------------------------------------------------------
# Makefile
#
# Simple makefile for building and installing land-surface-temperature.
#-----------------------------------------------------------------------------
.PHONY: all install clean

all:
	echo "make all in src..."; \
        (cd src; $(MAKE))

install: all
	echo "make install in scripts..."; \
        (cd scripts; $(MAKE) install)
	echo "make install in static_data..."; \
        (cd static_data; $(MAKE) install)
	echo "make install in src..."; \
        (cd src; $(MAKE) install)

clean:
	echo "make clean in src..."; \
        (cd src; $(MAKE) clean)

