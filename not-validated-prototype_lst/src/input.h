
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
    SAT_LANDSAT_4,
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
    char *scene_id;             /* SceneID */
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
    double bounding_coords[4];  /* Bounding Coordinates */
} Input_meta_t;


/* Structure for the 'input' data type */
typedef struct
{
    Input_meta_t meta;          /* Input metadata */
    int lines;
    int samples;
    float x_pixel_size;
    float y_pixel_size;
    char reference_band_name[30];
    char *band_name[MAX_INPUT_BANDS];
    FILE *band_fd[MAX_INPUT_BANDS];
    float scale_factor[MAX_INPUT_BANDS];
    int fill_value[MAX_INPUT_BANDS];
    float thermal_rad_gain;       /* Thermal radiance gain */
    float thermal_rad_bias;       /* Thermal radiance bias */
} Input_Data_t;


/* Prototypes */
Input_Data_t *open_input(Espa_internal_meta_t *metadata);

int close_input(Input_Data_t *input);

int read_input(Input_Data_t *input_data,
               float *band_thermal,
               int16_t *band_elevation,
               int pixel_count);

bool GetXMLInput(Input_Data_t *input,
                 Espa_internal_meta_t *metadata);


#endif
