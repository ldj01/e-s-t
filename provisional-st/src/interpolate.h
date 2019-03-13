#ifndef INTERPOLATE_H
#define INTERPOLATE_H

#include "grid_points.h"
#include "modtran_utils.h"


/* Defines index locations for the parameters in the at_height array */
typedef enum
{
    AHP_TRANSMISSION,
    AHP_UPWELLED_RADIANCE,
    AHP_DOWNWELLED_RADIANCE,
    AHP_NUM_PARAMETERS
} AT_HEIGHT_PARAMETERS;

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

void interpolate_to_height
(
    MODTRAN_POINT *modtran_points, /* I: results from MODTRAN runs */
    int *cell_vertices,            /* I: current cell vertices */
    double interpolate_to,         /* I: current landsat pixel height */
    double at_heights[][AHP_NUM_PARAMETERS]  /* I/O: interpolated parameters */
);

void interpolate_to_location
(
    GRID_POINTS *points,         /* I: The coordinate points */
    int *vertices,               /* I: The vertices for the points to use */
    double at_height[][AHP_NUM_PARAMETERS], /* I: current height atmospheric
                                                  results */
    double interpolate_easting,  /* I: interpolate to easting */
    double interpolate_northing, /* I: interpolate to northing */
    double *parameters           /* O: interpolated pixel atmospheric 
                                       parameters */
);

#endif
