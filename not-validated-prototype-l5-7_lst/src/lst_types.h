#ifndef LST_TYPES_H
#define LST_TYPES_H


#include <limits.h>


typedef struct
{
    char path[PATH_MAX];
    char command[PATH_MAX];
    float latitude;
    float longitude;
    float height;
} POINT_INFO;


typedef struct
{
    float ul_lat;
    float ul_lon;
    float lr_lat;
    float lr_lon;

    int min_eye;
    int max_eye;
    int min_jay;
    int max_jay;

    int num_eyes;
    int num_jays;
    int num_points;
    int num_modtran_runs;

    POINT_INFO *modtran_runs;

    float *eye;
    float *jay;
    float *lat;
    float *lon;
    float *utm_easting;
    float *utm_northing;
} REANALYSIS_POINTS;


#endif /* LST_TYPES_H */
