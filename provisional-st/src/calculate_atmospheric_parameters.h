#ifndef CALCULATE_ATMOSPHERIC_PARAMETERS_H
#define CALCULATE_ATMOSPHERIC_PARAMETERS_H

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


/* Function prototypes */

#endif /* CALCULATE_ATMOSPHERIC_PARAMETERS_H */

