#ifndef CALCULATE_POINT_ATMOSPHERIC_PARAMETERS_H
#define CALCULATE_POINT_ATMOSPHERIC_PARAMETERS_H


#include <stdbool.h>


#include "lst_types.h"
#include "input.h"


int calculate_point_atmospheric_parameters
(
    Input_t *input,        /* I: input structure */
    int num_points,        /* I: number of narr points */
    float albedo,          /* I: albedo */
    POINT_INFO *case_list, /* I: modtran run list */
    float **results,       /* O: atmospheric parameter for modtarn run */
    bool verbose           /* I: value to indicate if intermediate messages
                                 will be printed */
);


#endif /* CALCULATE_POINT_ATMOSPHERIC_PARAMETERS_H */
