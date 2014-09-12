#ifndef CFMASK_H
#define CFMASK_H
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>
#include "error.h"
#include "espa_metadata.h"
#include "parse_metadata.h"
#include "write_metadata.h"
#include "envi_header.h"
#include "espa_geoloc.h"
#include "raw_binary_io.h"

#define LST_VERSION "1.0.0"

typedef signed short int16;
typedef unsigned char uint8;

/* There are currently a maximum of 6 reflective bands in the output surface
   reflectance product */
#define NBAND_REFL_MAX 6

#define MAX_STR_LEN (510)

/* Input file type definition */
typedef enum {
  INPUT_TYPE_NULL = -1,
  INPUT_TYPE_BINARY = 0, 
  INPUT_TYPE_MAX
} Input_type_t;

/* Satellite type definition */
typedef enum {
  SAT_NULL = -1,
  SAT_LANDSAT_1 = 0, 
  SAT_LANDSAT_2, 
  SAT_LANDSAT_3, 
  SAT_LANDSAT_4, 
  SAT_LANDSAT_5, 
  SAT_LANDSAT_7, 
  SAT_LANDSAT_8, 
  SAT_MAX
} Sat_t;

/* Instrument type definition */
typedef enum {
  INST_NULL = -1,
  INST_MSS = 0, 
  INST_TM,
  INST_ETM, 
  INST_OLI_TIRS, 
  INST_MAX
} Inst_t;

/* World Reference System (WRS) type definition */
typedef enum {
  WRS_NULL = -1,
  WRS_1 = 0, 
  WRS_2,
  WRS_MAX
} Wrs_t;

/* Band gain settings (ETM+ only) */
typedef enum {
  GAIN_NULL = -1,
  GAIN_HIGH = 0, 
  GAIN_LOW, 
  GAIN_MAX
} Gain_t;

#endif
