#ifndef GRID_POINTS_H
#define GRID_POINTS_H

#include <stdint.h>


#define EQUATORIAL_RADIUS (6378137.0)

typedef struct
{
    int16_t index;
    int8_t run_modtran;
    uint8_t row;
    uint8_t col;
    uint16_t narr_row;
    uint16_t narr_col;
    float lon;          /* longitude (degrees) */
    float lat;          /* latitude (degrees) */
    float map_x;
    float map_y;
} GRID_POINT;

typedef struct
{
    int count;
    int rows;
    int cols;
    GRID_POINT *points;
} GRID_POINTS;

/* Defines the distance to the current pixel, along with the index of the
   point So that we can find the index of the closest point to start 
   determining the correct cell to use */
typedef struct
{
    int index;
    double distance;
} GRID_ITEM;


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
);

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
);

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
);

void free_grid_points
(
    GRID_POINTS *grid_points   /* I: Grid points to free */
);

int load_grid_points
(
    GRID_POINTS *grid_points
);

#endif
