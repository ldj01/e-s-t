#ifndef LST_TYPES_H
#define LST_TYPES_H


#include <stdbool.h>
#include <limits.h>


typedef struct
{
    char path[PATH_MAX];
    char command[PATH_MAX];
    double latitude;
    double longitude;
    double height;
    bool run_modtran;
} MODTRAN_INFO;


typedef struct
{
    double ul_lat;
    double ul_lon;
    double lr_lat;
    double lr_lon;

    int min_row;
    int max_row;
    int min_col;
    int max_col;

    int num_rows;
    int num_cols;
    int num_points;
    int num_modtran_runs;

    MODTRAN_INFO *modtran_runs;

    double *row;
    double *col;
    double *lat;
    double *lon;
    double *utm_easting;
    double *utm_northing;
    bool *use_point;
} REANALYSIS_POINTS;


#endif /* LST_TYPES_H */
