#ifndef LST_TYPES_H
#define LST_TYPES_H


#include <limits.h>


typedef struct
{
    char full_path[PATH_MAX];
    float latitude;
    float longitude;
    float height;
} POINT_INFO;


#endif /* LST_TYPES_H */
