#ifndef CONST_H
#define CONST_H


#include <math.h>


#define LST_VERSION "0.0.5"
#define LST_NO_DATA_VALUE (-9999.0)


#define LST_THERMAL_RADIANCE_PRODUCT_NAME "lst_thermal_radiance"
#define LST_THERMAL_RADIANCE_BAND_NAME "lst_thermal_radiance"
#define LST_THERMAL_RADIANCE_SHORT_NAME "LST_THERMAL_RADIANCE"
#define LST_THERMAL_RADIANCE_LONG_NAME "thermal band converted to radiance"
#define LST_RADIANCE_UNITS "radiance (W m^(-2) sr^(-1) mu^(-1))"

#define LST_ATMOS_TRANS_PRODUCT_NAME "lst_atmospheric_transmittance"
#define LST_ATMOS_TRANS_BAND_NAME "lst_atmospheric_transmittance"
#define LST_ATMOS_TRANS_SHORT_NAME "LST_ATMOSPHERIC_TRANSMITTANCE"
#define LST_ATMOS_TRANS_LONG_NAME "atmospheric transmittance"

#define LST_UPWELLED_RADIANCE_PRODUCT_NAME "lst_upwelled_radiance"
#define LST_UPWELLED_RADIANCE_BAND_NAME "lst_upwelled_radiance"
#define LST_UPWELLED_RADIANCE_SHORT_NAME "LST_UPWELLED_RADIANCE"
#define LST_UPWELLED_RADIANCE_LONG_NAME "upwelled radiance"

#define LST_DOWNWELLED_RADIANCE_PRODUCT_NAME "lst_downwelled_radiance"
#define LST_DOWNWELLED_RADIANCE_BAND_NAME "lst_downwelled_radiance"
#define LST_DOWNWELLED_RADIANCE_SHORT_NAME "LST_DOWNWELLED_RADIANCE"
#define LST_DOWNWELLED_RADIANCE_LONG_NAME "downwelled radiance"


#define TWO_PI (2.0 * PI)
#define HALF_PI (PI / 2.0)


#define DEGREES_PER_RADIAN (180.0 / PI)
#define RADIANS_PER_DEGREE (PI / 180.0)


#ifndef SUCCESS
    #define SUCCESS  0
#endif


#ifndef FAILURE
    #define FAILURE 1
#endif


#define NUM_ELEVATIONS 9
#define NARR_ROWS 277
#define NARR_COLS 349


#define MINSIGMA 1e-5


#define MAX_STR_LEN 512


#define UTM_SCALE_FACTOR (0.9996)
#define UTM_EQUATORIAL_RADIUS (6378137.0)
#define UTM_POLAR_RADIUS (6356752.3142)
#define UTM_FALSE_EASTING (500000.0)


/* Provides for the element size of the LST point result records, along with
   the locations of the individual elements.

   Only required because a 2d-array was used instead of a data structure.
*/
typedef enum {
    LST_LATITUDE = 0,
    LST_LONGITUDE,
    LST_HEIGHT,
    LST_TRANSMISSION,
    LST_UPWELLED_RADIANCE,
    LST_DOWNWELLED_RADIANCE,
    LST_NUM_ELEMENTS
} LST_RESULTS;


#endif
