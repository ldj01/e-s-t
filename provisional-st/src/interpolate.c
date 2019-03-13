/******************************************************************************
FILE: interpolate.c

PURPOSE:  This file includes routines to handle interpolation for the
          surface temperature atmospheric parameter calculations.
******************************************************************************/
#include <math.h>
#include "interpolate.h"

/******************************************************************************
METHOD:  interpolate_parameters

PURPOSE: Interpolate atmospheric parameters to the location of the
         current pixel.
******************************************************************************/
void interpolate_parameters
(
    MODTRAN_POINT *modtran_points, /* I: results from MODTRAN runs */
    GRID_POINTS *points,           /* I: coordinate points */
    int *cell_vertices,            /* I: current cell vertices */
    double interpolate_height,     /* I: current landsat pixel height */
    double interpolate_easting,    /* I: interpolate to easting */
    double interpolate_northing,   /* I: interpolate to northing */
    double *parameters             /* O: interpolated pixel atmospheric
                                         parameters */
)
{
    int parameter;
    int elevation;
    int below;
    int above;
    int vertex;
    double at_heights[NUM_CELL_POINTS][AHP_NUM_PARAMETERS];
    double below_parameters[AHP_NUM_PARAMETERS];
    double above_parameters[AHP_NUM_PARAMETERS];
    double w[NUM_CELL_POINTS];
    double total = 0.0;
    double inv_total;
    double interp_ratio; /* interpolation ratio */

    /* Interpolate three parameters to the height at each of the four
       closest points. */
    for (vertex = 0; vertex < NUM_CELL_POINTS; vertex++)
    {
        MODTRAN_POINT point = modtran_points[cell_vertices[vertex]];
        double *at_height = at_heights[vertex];

        /* Find the heights between which the pixel height sits.
           The heights are in increasing order. */
        for (elevation = 0; elevation < point.count; elevation++)
        {
            if (point.elevations[elevation].elevation >= interpolate_height)
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
            interp_ratio = (interpolate_height -
                            point.elevations[above].elevation)
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

    /* Interpolate parameters at appropriate height to location of the
       current pixel. */

    /* Shepard's method */
    for (vertex = 0; vertex < NUM_CELL_POINTS; vertex++)
    {
        double delta_x = points->points[cell_vertices[vertex]].map_x
                       - interpolate_easting;
        double delta_y = points->points[cell_vertices[vertex]].map_y
                       - interpolate_northing;

        w[vertex] = 1.0 / sqrt(delta_x*delta_x + delta_y*delta_y);

        total += w[vertex];
    }

    /* Normalize the weights for each vertex. */
    inv_total = 1/total;
    for (vertex = 0; vertex < NUM_CELL_POINTS; vertex++)
    {
        w[vertex] *= inv_total;
    }

    /* For each parameter apply each vertex's weighted value */
    for (parameter = 0; parameter < AHP_NUM_PARAMETERS; parameter++)
    {
        parameters[parameter] = 0.0;
        for (vertex = 0; vertex < NUM_CELL_POINTS; vertex++)
        {
            parameters[parameter] += w[vertex] * at_heights[vertex][parameter];
        }
    }
}
