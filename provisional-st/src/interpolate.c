/******************************************************************************
FILE: interpolate.c

PURPOSE:  This file includes routines to handle interpolation for the
          surface temperature atmospheric parameter calculations.
******************************************************************************/
#include <math.h>
#include "interpolate.h"

/******************************************************************************
METHOD:  interpolate_to_height

PURPOSE: Interpolate to height of current pixel
******************************************************************************/
void interpolate_to_height
(
    MODTRAN_POINT *modtran_points, /* I: results from MODTRAN runs */
    int *cell_vertices,            /* I: current cell vertices */
    double interpolate_to,         /* I: current landsat pixel height */
    double at_heights[][AHP_NUM_PARAMETERS]  /* I/O: interpolated parameters */
)
{
    int parameter;
    int elevation;
    int below;
    int above;

    double below_parameters[AHP_NUM_PARAMETERS];
    double above_parameters[AHP_NUM_PARAMETERS];

    double interp_ratio; /* interpolation ratio */

    int vertex;
    for (vertex = 0; vertex < NUM_CELL_POINTS; vertex++)
    {
        int current_index = cell_vertices[vertex];
        MODTRAN_POINT point = modtran_points[current_index];
        double *at_height = at_heights[vertex];

        /* Find the heights between which the pixel height sits.
           The heights are in increasing order. */
        for (elevation = 0; elevation < point.count; elevation++)
        {
            if (point.elevations[elevation].elevation >= interpolate_to)
                break;
        }
        if (elevation < point.count)
        {
            above = elevation;
            if (above != 0)
                below = above - 1;
            else
                /* All heights are above the pixel height, so use the first
                   value. */
                below = 0;
        }
        else
        {
            /* All heights are below the pixel height, so use the final
               value. */
            above = point.count - 1;
            below = above;
        }

        below_parameters[AHP_TRANSMISSION] =
            point.elevations[below].transmission;
        below_parameters[AHP_UPWELLED_RADIANCE] =
            point.elevations[below].upwelled_radiance;
        below_parameters[AHP_DOWNWELLED_RADIANCE] =
            point.elevations[below].downwelled_radiance;

        if (above == below)
        {
            /* Use the below parameters since the same */
            at_height[AHP_TRANSMISSION] =
                below_parameters[AHP_TRANSMISSION];
            at_height[AHP_UPWELLED_RADIANCE] =
                below_parameters[AHP_UPWELLED_RADIANCE];
            at_height[AHP_DOWNWELLED_RADIANCE] =
                below_parameters[AHP_DOWNWELLED_RADIANCE];
        }
        else
        {
            /* Interpolate between the heights for each parameter */
            interp_ratio = (interpolate_to - point.elevations[above].elevation)
                         / (point.elevations[above].elevation -
                            point.elevations[below].elevation);

            above_parameters[AHP_TRANSMISSION] =
                point.elevations[above].transmission;
            above_parameters[AHP_UPWELLED_RADIANCE] =
                point.elevations[above].upwelled_radiance;
            above_parameters[AHP_DOWNWELLED_RADIANCE] =
                point.elevations[above].downwelled_radiance;

            for (parameter = 0; parameter < AHP_NUM_PARAMETERS; parameter++)
            {
                at_height[parameter] = interp_ratio
                                     * (above_parameters[parameter] -
                                        below_parameters[parameter])
                                     + above_parameters[parameter];
            }
        }
    } /* vertex loop */
}


/******************************************************************************
METHOD:  interpolate_to_location

PURPOSE: Interpolate to location of current pixel
******************************************************************************/
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
)
{
    int point;
    int parameter;

    double w[NUM_CELL_POINTS];
    double total = 0.0;
    double inv_total;

    /* Shepard's method */
    for (point = 0; point < NUM_CELL_POINTS; point++)
    {
        double delta_x = points->points[vertices[point]].map_x
                       - interpolate_easting;
        double delta_y = points->points[vertices[point]].map_y
                       - interpolate_northing;

        w[point] = 1.0 / sqrt(delta_x*delta_x + delta_y*delta_y);

        total += w[point];
    }

    /* Normalize the weights for each vertex. */
    inv_total = 1/total;
    for (point = 0; point < NUM_CELL_POINTS; point++)
    {
        w[point] *= inv_total;
    }

    /* For each parameter apply each vertex's weighted value */
    for (parameter = 0; parameter < AHP_NUM_PARAMETERS; parameter++)
    {
        parameters[parameter] = 0.0;
        for (point = 0; point < NUM_CELL_POINTS; point++)
        {
            parameters[parameter] += (w[point] * at_height[point][parameter]);
        }
    }
}
