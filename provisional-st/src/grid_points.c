/*****************************************************************************
FILE: grid_points.c

PURPOSE:  This file contains routines to determine grid points.
*****************************************************************************/
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include "espa_geoloc.h"
#include "const.h"
#include "utilities.h"
#include "grid_points.h"


/******************************************************************************
METHOD:  qsort_grid_compare_function

PURPOSE: A qsort routine that can be used with the GRID_ITEM items to sort by
         distance

RETURN: int: -1 (a<b), 1 (b<a), 0 (a==b) 
******************************************************************************/
static int qsort_grid_compare_function
(
    const void *grid_item_a,
    const void *grid_item_b
)
{
    double a = (*(GRID_ITEM*)grid_item_a).distance;
    double b = (*(GRID_ITEM*)grid_item_b).distance;

    if (a < b)
        return -1;
    else if (b < a)
        return 1;

    return 0;
}


/******************************************************************************
METHOD:  haversine_distance

PURPOSE: Calculates the great-circle distance between 2 points in meters.
         The points are given in decimal degrees.  The Haversine formula
         is used. 

RETURN: double - The great-circle distance in meters between the points.

NOTE: This is based on the haversine_distance function in the ST Python
      scripts.
******************************************************************************/
static double haversine_distance
(
    double lon_1,  /* I: longitude (radians) for the first point */
    double lat_1,  /* I: latitude (radians) for the first point */
    double lon_2,  /* I: longitude (radians) for the second point */
    double lat_2   /* I: latitude (radians) for the second point */
)
{
    double sin_lon = sin((lon_2 - lon_1)*0.5);
    double sin_lat = sin((lat_2 - lat_1)*0.5);

    /* Compute and return the distance */
    return EQUATORIAL_RADIUS * 2
        + asin(sqrt(sin_lat*sin_lat + cos(lat_1)*cos(lat_2)*sin_lon*sin_lon));
}


/*****************************************************************************
METHOD:  determine_grid_point_distances

PURPOSE: Determines the distances for the current set of grid points.

NOTE: The indexes of the grid points are assumed to be populated.
*****************************************************************************/
void determine_grid_point_distances
(
    double *grid_pt_lon,       /* I: Grid point longitudes (radians) */
    double *grid_pt_lat,       /* I: Grid point latitudes (radians) */
    double longitude,          /* I: Longitude of the current line/sample
                                     (radians) */
    double latitude,           /* I: Latitude of the current line/sample
                                     (radians) */
    int num_grid_points,       /* I: The number of grid points to operate on */
    GRID_ITEM *grid_points     /* I/O: Sorted to determine the center grid
                                       point */
)
{
    int point;
    GRID_ITEM *pt = grid_points;  /* array pointer */

    /* Populate the distances to the grid points */
    for (point = 0; point < num_grid_points; point++, pt++)
    {
        pt->distance = haversine_distance(grid_pt_lon[pt->index],
                                          grid_pt_lat[pt->index],
                                          longitude, latitude);
    }
}


/*****************************************************************************
METHOD:  determine_center_grid_point

PURPOSE: Determines the index of the center point from the current set of grid
         points.

NOTE: The indexes of the grid points are assumed to be populated.

RETURN: type = int
    Value  Description
    -----  -------------------------------------------------------------------
    index  The index of the center point
*****************************************************************************/
int determine_center_grid_point
(
    double *grid_pt_lon,       /* I: Grid point longitudes (radians) */
    double *grid_pt_lat,       /* I: Grid point latitudes (radians) */
    double longitude,          /* I: Longitude of the current line/sample
                                     (radians) */
    double latitude,           /* I: Latitude of the current line/sample
                                     (radians) */
    int num_grid_points,       /* I: Number of grid points to operate on */
    GRID_ITEM *grid_points     /* I/O: Sorted to determine the center grid
                                       point */
)
{
    determine_grid_point_distances(grid_pt_lon, grid_pt_lat, longitude,
                                   latitude, num_grid_points, grid_points);

    /* Sort them to find the closest one */
    qsort (grid_points, num_grid_points, sizeof (GRID_ITEM),
           qsort_grid_compare_function);

    return grid_points[0].index;
}


/*****************************************************************************
METHOD:  determine_first_center_grid_point

PURPOSE: Determines the index of the first center point to use for the current
         line.  Only called when the fist valid point for a line is
         encountered.  The point is determined from all of the available
         points.

RETURN: type = int
    Value  Description
    -----  -------------------------------------------------------------------
    index  The index of the center point
*****************************************************************************/
int determine_first_center_grid_point
(
    double *grid_pt_lon,       /* I: Grid point longitudes (radians) */
    double *grid_pt_lat,       /* I: Grid point latitudes (radians) */
    double longitude,          /* I: Longitude of the current line/sample
                                     (radians) */
    double latitude,           /* I: Latitude of the current line/sample
                                     (radians) */
    int num_grid_points,       /* I: Number of grid points to operate on */
    GRID_ITEM *grid_points     /* I/O: Memory passed in, populated and
                                       sorted to determine the center grid
                                       point */
)
{
    int point;

    /* Assign the point indexes for all grid points */
    for (point = 0; point < num_grid_points; point++)
    {
        grid_points[point].index = point;
    }

    return determine_center_grid_point(grid_pt_lon, grid_pt_lat, longitude,
                                       latitude, num_grid_points, grid_points);
}


/*****************************************************************************
Method:  load_grid_points_hdr

Description:  Loads the grid points header information.

Notes:
    1. The grid point header file must be present in the current working 
       directory.

RETURN: SUCCESS
        FAILURE
*****************************************************************************/
int load_grid_points_hdr
(
    GRID_POINTS *grid_points
)
{
    char FUNC_NAME[] = "load_grid_points_hdr";
    FILE *grid_fd = NULL;
    int status;
    char header_filename[] = "grid_points.hdr";
    char errmsg[PATH_MAX];

    /* Open the grid header file. */
    grid_fd = fopen(header_filename, "r");
    if (grid_fd == NULL)
    {
        snprintf(errmsg, sizeof(errmsg), "Failed opening %s", header_filename);
        RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
    }

    /* Read the grid header file. */
    errno = 0;
    status = fscanf(grid_fd, "%d\n%d\n%d", &grid_points->count,
                    &grid_points->rows, &grid_points->cols);
    if (status != 3 || errno != 0)
    {
        fclose(grid_fd);
        snprintf(errmsg, sizeof(errmsg), "Failed reading %s", header_filename);
        RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
    }

    fclose(grid_fd);

    return SUCCESS;
}


/*****************************************************************************
Method:  load_grid_points

Description:  Loads the grid points into a data structure.

Notes:
    1. The grid point files must be present in the current working directory.

RETURN: SUCCESS
        FAILURE
*****************************************************************************/
int load_grid_points
(
    GRID_POINTS *grid_points
)
{
    char FUNC_NAME[] = "load_grid_points";

    FILE *grid_fd = NULL;

    int status;

    char binary_filename[] = "grid_points.bin";
    char errmsg[PATH_MAX];

    /* Initialize the points */
    grid_points->points = NULL;

    if (load_grid_points_hdr(grid_points) != SUCCESS)
    {
        RETURN_ERROR("Failed loading grid point header information",
                     FUNC_NAME, FAILURE);
    }

    grid_points->points = malloc(grid_points->count * sizeof(GRID_POINT));
    if (grid_points->points == NULL)
    {
        RETURN_ERROR("Failed allocating memory for grid points",
                     FUNC_NAME, FAILURE);
    }

    /* Open the grid point file */
    grid_fd = fopen(binary_filename, "rb");
    if (grid_fd == NULL)
    {
        snprintf(errmsg, sizeof(errmsg), "Failed opening %s", binary_filename);
        RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
    }

    /* Read the grid points */
    status = fread(grid_points->points, sizeof(GRID_POINT),
                   grid_points->count, grid_fd);
    if (status != grid_points->count || errno != 0)
    {
        fclose(grid_fd);
        snprintf(errmsg, sizeof(errmsg), "Failed reading %s", binary_filename);
        RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
    }

    fclose(grid_fd);

    return SUCCESS;
}


/*****************************************************************************
Method:  free_grid_points

Description:  Free allocated memory for the grid points.
*****************************************************************************/
void free_grid_points
(
    GRID_POINTS *grid_points
)
{
    free(grid_points->points);
    grid_points->points = NULL;
}
