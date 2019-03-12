/*****************************************************************************
FILE: modtran_utils.c

PURPOSE:  This file contains routines for working with MODTRAN data.
*****************************************************************************/
#include <stdio.h>
#include "espa_geoloc.h"
#include "const.h"
#include "utilities.h"

#include "modtran_utils.h"

/*****************************************************************************
Method:  load_elevations

Description:  Loads the grid elevations into a data structure.

Notes:
    1. The grid elevation file must be present in the current working directory.
    2. The grid elevation entries should be in sync with the grid file.

RETURN: SUCCESS
        FAILURE
*****************************************************************************/
static int load_elevations
(
    MODTRAN_POINTS *modtran_points
)
{
    char FUNC_NAME[] = "load_elevations";

    FILE *elevation_fd = NULL;

    int status;
    int index;   /* Index into point structure */

    char elevation_filename[] = "grid_elevations.txt";
    char errmsg[PATH_MAX];

    snprintf(errmsg, sizeof(errmsg), "Failed reading %s", elevation_filename);

    elevation_fd = fopen(elevation_filename, "r");
    if (elevation_fd == NULL)
    {
        RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
    }

    /* Read the elevations into the 0 elevation positions in the MODTRAN 
       point structure.  The file and structure should have the same order. */
    for (index = 0; index < modtran_points->count; index++)
    {
        MODTRAN_POINT *modtran_ptr = &modtran_points->points[index];

        /* Keep looking for a modtran point that was actually run. */
        if (modtran_ptr->ran_modtran == 0)
        {
            continue; 
        }

        status = fscanf(elevation_fd, "%lf %lf\n", 
            &(modtran_ptr->elevations->elevation),
            &(modtran_ptr->elevations->elevation_directory));
        if (status <= 0)
        {
            RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
        }
    }

    fclose(elevation_fd);

    return SUCCESS;
}


/*****************************************************************************
Method:  initialize_modtran_points

Description:  Allocate the memory need to hold the MODTRAN results and
              initialize known values.

RETURN: SUCCESS
        FAILURE
*****************************************************************************/
int initialize_modtran_points
(
    GRID_POINTS *grid_points,      /* I: The coordinate points */
    MODTRAN_POINTS *modtran_points /* O: Memory Allocated */
)
{
    char FUNC_NAME[] = "initialize_modtran_points";

    double gndalt[MAX_NUM_ELEVATIONS];

    int index;
    int num_elevations;            /* Number of elevations actually used */
    int elevation_index;           /* Index into elevations */
    int status;                    /* Function return status */

    FILE *modtran_elevation_fd = NULL;

    char modtran_elevation_filename[] = "modtran_elevations.txt";
    char errmsg[PATH_MAX];

    snprintf(errmsg, sizeof(errmsg), "Failed reading %s", 
        modtran_elevation_filename);

    modtran_elevation_fd = fopen(modtran_elevation_filename, "r");
    if (modtran_elevation_fd == NULL)
    {
        RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
    }

    status = fscanf(modtran_elevation_fd, "%d\n", &num_elevations);
    if (status <= 0)
    {
        RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
    }

    /* Read the elevations into the gndalt structure. */ 
    for (index = 0; index < num_elevations; index++)
    {
        status = fscanf(modtran_elevation_fd, "%lf\n", &gndalt[index]);
        if (status <= 0)
        {
            RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
        }
    }

    fclose(modtran_elevation_fd);

    modtran_points->count = grid_points->count;

    modtran_points->points = malloc(modtran_points->count *
                                    sizeof(MODTRAN_POINT));
    if (modtran_points->points == NULL)
    {
        RETURN_ERROR("Failed allocating memory for modtran points",
                     FUNC_NAME, FAILURE);
    }

    for (index = 0; index < modtran_points->count; index++)
    {
        /* convenience pointers */
        MODTRAN_POINT *modtran_ptr = &modtran_points->points[index];
        GRID_POINT *grid_ptr = &grid_points->points[index];
        
        modtran_ptr->count = num_elevations;
        modtran_ptr->ran_modtran = grid_ptr->run_modtran;
        modtran_ptr->row = grid_ptr->row;
        modtran_ptr->col = grid_ptr->col;
        modtran_ptr->narr_row = grid_ptr->narr_row;
        modtran_ptr->narr_col = grid_ptr->narr_col;
        modtran_ptr->lon = grid_ptr->lon;
        modtran_ptr->lat = grid_ptr->lat;
        modtran_ptr->map_x = grid_ptr->map_x;
        modtran_ptr->map_y = grid_ptr->map_y;

        modtran_ptr->elevations = malloc(num_elevations*
                                         sizeof(MODTRAN_ELEVATION));
        if (modtran_ptr->elevations == NULL)
        {
            RETURN_ERROR("Failed allocating memory for modtran point"
                         " elevations", FUNC_NAME, FAILURE);
        }

        /* Iterate over the elevations and assign the elevation values. */
        for (elevation_index = 0; elevation_index < num_elevations; 
            elevation_index++)
        {
            MODTRAN_ELEVATION *mt_elev =
                                    &modtran_ptr->elevations[elevation_index];

            mt_elev->elevation = gndalt[elevation_index];
            mt_elev->elevation_directory = gndalt[elevation_index];
        }
    }

    /* Load the first elevation values if needed. */
    if (load_elevations(modtran_points) != SUCCESS)
    {
        RETURN_ERROR("calling load_elevations", FUNC_NAME, EXIT_FAILURE);
    }

    return SUCCESS;
}


/*****************************************************************************
Method:  free_modtran_points

Description:  Free allocated memory for the MODTRAN points.
*****************************************************************************/
void free_modtran_points
(
    MODTRAN_POINTS *modtran_points
)
{
    int index;     /* Index into MODTRAN points structure */

    for (index = 0; index < modtran_points->count; index++)
    {
        free(modtran_points->points[index].elevations);
        modtran_points->points[index].elevations = NULL;
    }

    free(modtran_points->points);
    modtran_points->points = NULL;
}
