
#ifndef BUILD_POINTS_H
#define BUILD_POINTS_H


#include "lst_types.h"
#include "input.h"


int build_points
(
    Input_t *input,           /* I: input structure */
    REANALYSIS_POINTS *points /* O: The coordinate points to be used */
);


void free_points_memory
(
    REANALYSIS_POINTS *points /* I: The coordinate points */
);


#endif /* BUILD_POINTS_H */
