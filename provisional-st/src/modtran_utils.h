#ifndef MODTRAN_UTILS_H
#define MODTRAN_UTILS_H

#include "grid_points.h"

#define MAX_NUM_ELEVATIONS 9

typedef struct
{
    double elevation;
    double elevation_directory;
    double transmission;
    double upwelled_radiance;
    double downwelled_radiance;
} MODTRAN_ELEVATION;


typedef struct
{
    int count;
    int ran_modtran;
    uint8_t row;
    uint8_t col;
    uint16_t narr_row;
    uint16_t narr_col;
    double lon;                    /* longitude (degrees) */
    double lat;                    /* latitude (degrees) */
    double map_x;
    double map_y;
    MODTRAN_ELEVATION *elevations;
} MODTRAN_POINT;


typedef struct
{
    int count;
    MODTRAN_POINT *points;
} MODTRAN_POINTS;


int initialize_modtran_points
(
    GRID_POINTS *grid_points,      /* I: The coordinate points */
    MODTRAN_POINTS *modtran_points /* O: Memory Allocated */
);

void free_modtran_points
(
    MODTRAN_POINTS *modtran_points
);

#endif
