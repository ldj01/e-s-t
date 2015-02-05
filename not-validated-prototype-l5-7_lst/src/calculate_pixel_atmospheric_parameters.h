#ifndef LST_H
#define LST_H


#include <stdbool.h>


#include "lst_types.h"
#include "input.h"


int third_pixels_post
(
    Input_t *input,          /* I: input structure */
    int num_points,          /* I: number of narr points */
    char *dem_infile,        /* I: address of input DEM filename */
    char *emi_infile,        /* I: address of input Emissivity filename */
    float **modtran_results, /* I: atmospheric parameter for modtarn run */
    bool verbose             /* I: value to indicate if intermediate messages
                                   will be printed */
);


#endif /* LST_H */
