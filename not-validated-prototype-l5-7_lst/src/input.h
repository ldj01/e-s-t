#ifndef INPUT_H
#define INPUT_H


#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <ctype.h>


#include "date.h"
#include "espa_metadata.h"
#include "parse_metadata.h"
#include "write_metadata.h"
#include "envi_header.h"
#include "espa_geoloc.h"
#include "raw_binary_io.h"


#define INPUT_FILL (0)
#define ANGLE_FILL (-999.0)
#define WRS_FILL (-1)
#define GAIN_BIAS_FILL (-999.0)


/* Satellite type definition */
typedef enum
{
    SAT_NULL = -1,
    SAT_LANDSAT_5,
    SAT_LANDSAT_7,
    SAT_LANDSAT_8,
    SAT_MAX
} Satellite_t;


/* Instrument type definition */
typedef enum
{
    INST_NULL = -1,
    INST_TM,
    INST_ETM,
    INST_OLI_TIRS,
    INST_MAX
} Instrument_t;


/* Band gain settings (ETM+ only) */
typedef enum
{
    GAIN_NULL = -1,
    GAIN_HIGH = 0,
    GAIN_LOW,
    GAIN_MAX
} Gain_t;


/* Structure for the metadata */
typedef struct
{
    Satellite_t satellite;      /* Satellite */
    Instrument_t instrument;    /* Instrument */
    Date_t acq_date;            /* Acq. date/time (scene center) */
    int zone;                   /* UTM zone number */
    Map_coord_t ul_map_corner;  /* Map projection coordinates of the upper
                                   left corner of the pixel in the upper left
                                   corner of the image */
    Map_coord_t lr_map_corner;  /* Map projection coordinates of the lower
                                   right corner of the pixel in the lower right
                                   corner of the image */
    Geo_coord_t ul_geo_corner;  /* Geo projection coordinates of the upper
                                   left corner of the pixel in the upper left
                                   corner of the image */
    Geo_coord_t lr_geo_corner;  /* Geo projection coordinates of the lower
                                   right corner of the pixel in the lower right
                                   corner of the image */
} Input_meta_t;


/* Structure for the thermal band data */
typedef struct
{
    char *filename;       /* Name of the image file */
    FILE *fd;             /* File pointer for file */
    bool is_open;         /* Has the file been opened flag */
    int band_index;       /* Index in the metadata */
    Img_coord_int_t size; /* (line/sample) size */
    float toa_gain;       /* TOA gain */
    float toa_bias;       /* TOA bias */
    int fill_value;       /* Fill value */
    float pixel_size[2];  /* Pixel size (x,y) */
} Input_band_t;


/* Structure for the 'input' data type */
typedef struct
{
    Input_meta_t meta;          /* Input metadata */
    Input_band_t thermal;       /* Input thermal band data */
} Input_t;


/* Prototypes */
Input_t *OpenInput (Espa_internal_meta_t *metadata);


bool GetInputThermLine (Input_t *input, int iline, float *thermal_data);


bool CloseInput (Input_t *input);


bool FreeInput (Input_t *input);


bool GetXMLInput (Input_t *input, Espa_internal_meta_t * metadata);


void split_filename
(
    const char *filename, /* I: Name of file to split */
    char *directory,      /* O: Directory portion of file name */
    char *scene_name,     /* O: Scene name portion of the file name */
    char *extension       /* O: Extension portion of the file name */
);


#endif
