
#ifndef CALCULATE_ATMOSPHERIC_PARAMETERS_H
#define CALCULATE_ATMOSPHERIC_PARAMETERS_H

#include "st_types.h"
#include "input.h"
#include "espa_geoloc.h"

#define L4_TM_SRS_COUNT (171)
#define L5_TM_SRS_COUNT (171)
#define L7_TM_SRS_COUNT (125)
#define L8_OLITIRS_SRS_COUNT (101)
#define MAX_SRS_COUNT (L5_TM_SRS_COUNT)
/* This emissivity/albedo is for water */
#define WATER_ALBEDO (0.1)
#define WATER_EMISSIVITY (1.0 - WATER_ALBEDO)
#define INV_WATER_ALBEDO (1.0 / WATER_ALBEDO)

#define INV_TWO (0.5)
#define INV_SIX (1.0 / 6.0)

/* Defines the index for the intermediate bands which are generated for the
 *    follow-on processing */
typedef enum
{
    BAND_TRANSMISSION,
    BAND_UPWELLED_RADIANCE,
    BAND_DOWNWELLED_RADIANCE,
    BAND_THERMAL,
    NUM_INTERMEDIATE_DATA_BANDS
} INTERMEDIATE_DATA_BANDS;


/* Defines the distance to the current pixel, along with the index of the
   point So that we can find the index of the closest point to start 
   determining the correct cell to use */
typedef struct
{
    int index;
    double distance;
} GRID_ITEM;


/* Defines index locations in the vertices array for the current cell to be
   used for interpolation of the pixel */
typedef enum
{
    LL_POINT,
    UL_POINT,
    UR_POINT,
    LR_POINT,
    NUM_CELL_POINTS
} CELL_POINTS;


/* Defines index locations for the parameters in the at_height array */
typedef enum
{
    AHP_TRANSMISSION,
    AHP_UPWELLED_RADIANCE,
    AHP_DOWNWELLED_RADIANCE,
    AHP_NUM_PARAMETERS
} AT_HEIGHT_PARAMETERS;


/* Function prototypes */

int calculate_point_atmospheric_parameters
(
    Input_Data_t *input,       /* I: Input structure */
    GRID_POINTS *grid_points,  /* I: The coordinate points */
    MODTRAN_POINTS *modtran_results /* I/O: Atmospheric parameters from
                                   modtran */
);

int calculate_pixel_atmospheric_parameters
(
    Input_Data_t *input,       /* I: input structure */
    GRID_POINTS *points,       /* I: The coordinate points */
    char *xml_filename,        /* I: XML filename */
    Espa_internal_meta_t xml_metadata, /* I: XML metadata */
    MODTRAN_POINTS *modtran_results /* I: results from MODTRAN runs */
);

void free_grid_points
(
    GRID_POINTS *grid_points   /* I: Grid points to free */
);

#endif /* CALCULATE_ATMOSPHERIC_PARAMETERS_H */

