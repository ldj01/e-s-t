
#ifndef CONST_H
#define CONST_H


#include <math.h>


#define LST_VERSION "0.0.4"
#define LST_NO_DATA_VALUE (-9999.0)


#define LST_THERMAL_RADIANCE_PRODUCT_NAME "lst_temp"
#define LST_THERMAL_RADIANCE_BAND_NAME "lst_thermal_radiance"
#define LST_THERMAL_RADIANCE_SHORT_NAME "LST_THERMAL_RADIANCE"
#define LST_THERMAL_RADIANCE_LONG_NAME "thermal band converted to radiance"
#define LST_RADIANCE_UNITS "radiance (W m^(-2) sr^(-1) mu^(-1))"

#define LST_ATMOS_TRANS_PRODUCT_NAME "lst_temp"
#define LST_ATMOS_TRANS_BAND_NAME "lst_atmospheric_transmittance"
#define LST_ATMOS_TRANS_SHORT_NAME "LST_ATMOSPHERIC_TRANSMITTANCE"
#define LST_ATMOS_TRANS_LONG_NAME "atmospheric transmittance"

#define LST_UPWELLED_RADIANCE_PRODUCT_NAME "lst_temp"
#define LST_UPWELLED_RADIANCE_BAND_NAME "lst_upwelled_radiance"
#define LST_UPWELLED_RADIANCE_SHORT_NAME "LST_UPWELLED_RADIANCE"
#define LST_UPWELLED_RADIANCE_LONG_NAME "upwelled radiance"

#define LST_DOWNWELLED_RADIANCE_PRODUCT_NAME "lst_temp"
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


#define EQUATORIAL_RADIUS (6378137.0)


/* Provides for the element size of the LST point result records, along with
   the locations of the individual elements.

   Only required because a 2d-array was used instead of a data structure.
*/
typedef enum
{
    MGPE_LATITUDE = 0,
    MGPE_LONGITUDE,
    MGPE_HEIGHT,
    MGPE_TRANSMISSION,
    MGPE_UPWELLED_RADIANCE,
    MGPE_DOWNWELLED_RADIANCE,
    MGPE_NUM_ELEMENTS
} MODTRAN_GRID_POINT_ELEMENTS;


/* Defines index locations in the array of grid points that specify the nine
   closest grid points */
typedef enum
{
    CC_GRID_POINT,
    LL_GRID_POINT,
    LC_GRID_POINT,
    UL_GRID_POINT,
    UC_GRID_POINT,
    UR_GRID_POINT,
    RC_GRID_POINT,
    LR_GRID_POINT,
    DC_GRID_POINT,
    NUM_GRID_POINTS
} GRID_POINT_INDEXES;


/* These are used in arrays, and they are position dependent */
typedef enum
{
    I_BAND_THERMAL,
    I_BAND_ELEVATION, /* This band and above are all from the XML */
    MAX_INPUT_BANDS
} Input_Bands_e;


#endif /* CONST_H */
