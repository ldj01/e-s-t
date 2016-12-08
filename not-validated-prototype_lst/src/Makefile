#-----------------------------------------------------------------------------
# Makefile
#
# For building land-surface-temperature.
#-----------------------------------------------------------------------------
.PHONY: all install clean

# Inherit from upper-level make.config
TOP = ../..
include $(TOP)/make.config

# Set up compile options
CC = gcc
RM = rm
EXTRA = -Wall $(EXTRA_OPTIONS)

# Define the include files
INC1 = utilities.h 2d_array.h calculate_atmospheric_parameters.h input.h output.h intermediate_data.h
INCDIR  = -I. -I$(XML2INC) -I$(ESPAINC)
NCFLAGS = $(EXTRA) $(INCDIR)

# Define the source code and object files
SRC1 = \
      utilities.o                              \
      2d_array.c                               \
      input.c                                  \
      output.c                                 \
      intermediate_data.c                      \
      calculate_atmospheric_parameters.c
OBJ1 = $(SRC1:.c=.o)

# Define the object libraries
EXLIB = -L$(ESPALIB) -l_espa_raw_binary -l_espa_common \
        -L$(XML2LIB) -lxml2 \
        -L$(LZMALIB) -llzma \
        -L$(ZLIBLIB) -lz
MATHLIB = -lm
LOADLIB = $(EXLIB) $(MATHLIB)

# Define the executable
EXE1 = lst_atmospheric_parameters

# Target for the executable
all: $(EXE1)

$(EXE1): $(OBJ1) $(INC1)
	$(CC) $(EXTRA) -o $(EXE1) $(OBJ1) $(LOADLIB)

install:
	install -d $(link_path)
	install -d $(lst_install_path)
	install -m 755 $(EXE1) $(lst_install_path)
	ln -sf $(lst_link_source_path)/$(EXE1) $(link_path)/$(EXE1)

clean:
	$(RM) -f *.o $(EXE1)

$(OBJ1): $(INC1)

.c.o:
	$(CC) $(NCFLAGS) -c $<

