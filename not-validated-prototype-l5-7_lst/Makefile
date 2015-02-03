#------------------------------------------------------------------------------
# Makefile
#
# Simple makefile for building and installing land-surface-temperature.
#------------------------------------------------------------------------------

SUBDIRS	= scripts src static_data

all:
	@for dir in $(SUBDIRS); do \
        echo "make all in $$dir..."; \
        (cd $$dir; $(MAKE) -f Makefile); done

install: all
	@for dir in $(SUBDIRS); do \
        echo "make install in $$dir..."; \
        (cd $$dir; $(MAKE) -f Makefile install); done

clean:
	@for dir in $(SUBDIRS); do \
        echo "make clean in $$dir..."; \
        (cd $$dir; $(MAKE) -f Makefile clean); done

