
#ifndef LST_TYPES_H
#define LST_TYPES_H

#include <limits.h>

typedef struct {
    int8_t index;
    int8_t run_modtran;
    uint8_t row;
    uint8_t col;
    uint16_t narr_row;
    uint16_t narr_col;
    float lon;
    float lat;
    float map_x;
    float map_y;
} GRID_POINT;


typedef struct {
    int count;
    int rows;
    int cols;
    GRID_POINT *points;
} GRID_POINTS;


typedef struct {
    double elevation;
    double transmission;
    double upwelled_radiance;
    double downwelled_radiance;
} MODTRAN_ELEVATION;


typedef struct {
    int count;
    int ran_modtran;
    uint8_t row;
    uint8_t col;
    uint16_t narr_row;
    uint16_t narr_col;
    double lon;
    double lat;
    double map_x;
    double map_y;
    MODTRAN_ELEVATION *elevations;
} MODTRAN_POINT;


typedef struct {
    int count;
    MODTRAN_POINT *points;
} MODTRAN_POINTS;

#endif /* LST_TYPES_H */
