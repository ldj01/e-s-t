#ifndef INPUT_H
#define INPUT_H


#include <stdlib.h>
#include <stdio.h>
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
#define NBAND_THM_MAX 2
#define LST_VERSION "1.0.0"
#define NUM_ELEVATIONS 9
#define NARR_ROW 349
#define NARR_COL 277
#define P_LAYER 29
#define STAN_LAYER 30
#define MAX_MODTRAN_LAYER 150


typedef signed short int16;
typedef unsigned char uint8;


typedef enum
{
    HDF_FILE = 0,
    BINARY_FILE
} File_type;


/* Input file type definition */
typedef enum
{
    INPUT_TYPE_NULL = -1,
    INPUT_TYPE_BINARY = 0,
    INPUT_TYPE_MAX
} Input_type_t;


/* Satellite type definition */
typedef enum
{
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
typedef enum
{
    INST_NULL = -1,
    INST_MSS = 0,
    INST_TM,
    INST_ETM,
    INST_OLI_TIRS,
    INST_MAX
} Inst_t;


/* World Reference System (WRS) type definition */
typedef enum
{
    WRS_NULL = -1,
    WRS_1 = 0,
    WRS_2,
    WRS_MAX
} Wrs_t;


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
    Sat_t sat;                  /* Satellite */
    Inst_t inst;                /* Instrument */
    Date_t acq_date;            /* Acq. date/time (scene center) */
    bool time_fill;             /* Acq. time fill; true = fill value (0h) */
    Date_t prod_date;           /* Prod date (must be available for ETM) */
    float sun_zen;              /* Solar zenith angle (deg; scene center) */
    float sun_az;               /* Solar azimuth angle (deg; scene center) */
    Wrs_t wrs_sys;              /* WRS system */
    int ipath;                  /* WRS path number */
    int irow;                   /* WRS row number */
    unsigned char fill;         /* Fill value */
    int iband_th[NBAND_THM_MAX];/* thermal band numbers */
    int therm_satu_value_ref;   /* saturation value of thermal product */
    int therm_satu_value_max;   /* maximum bt value (degrees Celsius) */
    Gain_t gain_setting_th;     /* Band gain settings Thermal */
    float scale_factor_th;      /* Scale factor for the thermal band */
    int zone;                   /* UTM zone number */
    float pixel_size[2];        /* pixel size (x,y) */
    float gain_th;              /* Thermal band gain */
    float bias_th;              /* Thermal band bias */
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


/* Structure for the 'input' data type */
typedef struct
{
    Input_type_t file_type;     /* Type of the input image files */
    Input_meta_t meta;          /* Input metadata */
    int nband_th;               /* Num of thermal input image files (0 or 1) */
    Img_coord_int_t size_th;    /* Input (thermal) file size */
    char *file_name_th;         /* Name of the thermal input image files */
    bool open_th;               /* thermal open flag */
    FILE *fp_bin_th;            /* File pointer for thermal binary file */
    int16 *therm_buf;           /* Input data buf (one line of therm data) */
} Input_t;


/* Prototypes */
Input_t *OpenInput (Espa_internal_meta_t *metadata);


bool GetInputThermLine (Input_t *this, int iline);


bool CloseInput (Input_t *this);


bool FreeInput (Input_t *this);


bool GetXMLInput (Input_t *this, Espa_internal_meta_t * metadata);


void split_filename
(
    const char *filename, /* I: Name of file to split */
    char *directory,      /* O: Directory portion of file name */
    char *scene_name,     /* O: Scene name portion of the file name */
    char *extension       /* O: Extension portion of the file name */
);


int get_args
(
    int argc,                 /* I: number of cmd-line args */
    char *argv[],             /* I: string of cmd-line args */
    char **xml_infile,        /* I: address of input XML metadata filename  */
    char **dem_infile,        /* I: address of input DEM filename */
    char **emissivity_infile, /* I: address of input emissivity filename */
    bool *verbose             /* O: verbose flag */
);


void usage ();


#endif
