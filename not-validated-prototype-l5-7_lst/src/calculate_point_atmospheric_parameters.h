#ifndef CALCULATE_POINT_ATMOSPHERIC_PARAMETERS_H
#define CALCULATE_POINT_ATMOSPHERIC_PARAMETERS_H


#include <stdbool.h>


#include "lst_types.h"
#include "input.h"


int calculate_point_atmospheric_parameters
(
    Input_t *input,            /* I: input structure */
    REANALYSIS_POINTS *points, /* I: The coordinate points */
    float albedo,              /* I: albedo */
    float **results,           /* O: atmospheric parameter for modtarn run */
    bool verbose               /* I: value to indicate if intermediate
                                     messages will be printed */
);


#endif /* CALCULATE_POINT_ATMOSPHERIC_PARAMETERS_H */
