#ifndef CALCULATE_PIXEL_ATMOSPHERIC_PARAMETERS_H
#define CALCULATE_PIXEL_ATMOSPHERIC_PARAMETERS_H


#include <stdbool.h>


#include "lst_types.h"
#include "input.h"


int calculate_pixel_atmospheric_parameters
(
    Input_t *input,            /* I: input structure */
    REANALYSIS_POINTS *points, /* I: The coordinate points */
    char *dem_infile,          /* I: address of input DEM filename */
    char *emi_infile,          /* I: address of input Emissivity filename */
    float **modtran_results,   /* I: atmospheric parameter for modtarn run */
    bool verbose               /* I: value to indicate if intermediate
                                     messages will be printed */
);


#endif /* CALCULATE_PIXEL_ATMOSPHERIC_PARAMETERS_H */
