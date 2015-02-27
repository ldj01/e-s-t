#ifndef CALCULATE_PIXEL_ATMOSPHERIC_PARAMETERS_H
#define CALCULATE_PIXEL_ATMOSPHERIC_PARAMETERS_H


#include <stdbool.h>


#include "lst_types.h"
#include "input.h"


int calculate_pixel_atmospheric_parameters
(
    Input_t *input,            /* I: input structure */
    REANALYSIS_POINTS *points, /* I: The coordinate points */
    char *xml_filename,        /* I: XML filename */
    char *dem_filename,        /* I: DEM filename */
    char *emi_filename,        /* I: Emissivity filename */
    float **modtran_results,   /* I: atmospheric parameter for MODTRAN run */
    bool verbose               /* I: value to indicate if intermediate
                                     messages will be printed */
);


#endif /* CALCULATE_PIXEL_ATMOSPHERIC_PARAMETERS_H */
