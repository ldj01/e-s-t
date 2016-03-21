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
INC = const.h utilities.h date.h 2d_array.h get_args.h input.h output.h \
      build_points.h build_modtran_input.h \
      calculate_point_atmospheric_parameters.h \
      calculate_pixel_atmospheric_parameters.h
INCDIR  = -I. -I$(XML2INC) -I$(ESPAINC)
NCFLAGS = $(EXTRA) $(INCDIR)

# Define the source code and object files
SRC = \
      utilities.c                              \
      2d_array.c                               \
      date.c                                   \
      input.c                                  \
      intermediate_data.c                      \
      output.c                                 \
      get_args.c                               \
      build_points.c                           \
      build_modtran_input.c                    \
      calculate_point_atmospheric_parameters.c \
      calculate_pixel_atmospheric_parameters.c \
      lst.c
OBJ = $(SRC:.c=.o)

# Define the object libraries
EXLIB = -L$(ESPALIB) -l_espa_raw_binary -l_espa_common \
        -L$(XML2LIB) -lxml2 \
        -L$(LZMALIB) -llzma \
        -L$(ZLIBLIB) -lz
MATHLIB = -lm
LOADLIB = $(EXLIB) $(MATHLIB)

# Define the executable
EXE = lst_intermediate_data

# Target for the executable
all: $(EXE)

$(EXE): $(OBJ) $(INC)
	$(CC) $(EXTRA) -o $(EXE) $(OBJ) $(LOADLIB)

install:
	install -d $(link_path)
	install -d $(lst_install_path)
	install -m 755 $(EXE) $(lst_install_path)
	ln -sf $(lst_link_source_path)/$(EXE) $(link_path)/$(EXE)

clean:
	$(RM) -f *.o $(EXE)

$(OBJ): $(INC)

.c.o:
	$(CC) $(NCFLAGS) -c $<

