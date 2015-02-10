
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>


#include "const.h"
#include "2d_array.h"
#include "utilities.h"
#include "input.h"
#include "lst_types.h"
#include "build_points.h"


int qsort_float_compare_function
(
    const void *a,
    const void *b
)
{
    if (*(float*)a < *(float*)b)
        return -1;
    else if (*(float*)b < *(float*)a)
        return 1;

    return 0;
}


/******************************************************************************
MODULE:  distance_in_utm

PURPOSE: Calculate distances between UTM coordiantes

RETURN: NONE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
10/10/2014   Song Guo         Original Development
******************************************************************************/
void distance_in_utm
(
    float e1,
    float n1,
    float e2,
    float n2,
    float *distance
)
{
    float s = 0.9996;           /* scale factor */
    float r = 6378137;          /* earth radius */
    float sr1, sr2, sr3;
    float edist;

    sr1 = s / (cos (e1 / r));
    sr2 = s / (cos (((e2 - e1) / 6.0) / r));
    sr3 = s / (cos (e2 / r));

    edist = ((e2 - e1) / 6.0) * (sr1 + 4 * sr2 + sr3);

    *distance = sqrt (edist * edist + (n2 - n1) * (n2 - n1));
}


/******************************************************************************
MODULE:  convert_ll_utm

PURPOSE: convert digital counts to radiance for thermal band
         [unit: W m^(-2) sr^(-1) mu^(-1) ]

RETURN: NONE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/30/2014   Song Guo         Original Development
******************************************************************************/
void convert_ll_utm
(
    Input_t *input,
    float *lat,
    float *lon,
    int num_points,
    float **narr_utm
)
{
    int i;
    float a = 6378137.0;        /* equatorial radius */
    float b = 6356752.3;        /* polar radius */
    float k0 = 0.9996;          /* scale factor */
    float e;                    /* eccentricity */
    float eprimesqrd;
    float n;
//    float rho;
    float nu;
    float a0;
    float b0;
    float c0;
    float d0;
    float e0;
    float ki;
    float kii;
    float kiii;
    float kiv;
    float kv;
    float zone_cm;
//    float zone_cm_rad;
    float delta_lon;
    float p;
    float lat_rad;
//    float lon_rad;
    float s;
    float northing;
    float easting;

    /* calculate zone central meridian in degrees and radians */
    zone_cm = (float) (6 * input->meta.zone - 183);
//    zone_cm_rad =  zone_cm * PI / 180.0;

    e = sqrt (1.0 - pow (b / a, 2));
    eprimesqrd = e * e / (1.0 - e * e);
    n = (a - b) / (a + b);

    /* calculate meridional arc length */
    a0 = a * (1.0 - n + (5.0 * n * n / 4.0) * (1.0 - n) +
              (81.0 * (float) pow (n, 4) / 64.0) * (1.0 - n));
    b0 = (3.0 * a * n / 2.0) * (1.0 - n - (7.0 * n * n / 8.0) * (1.0 - n) +
                                55.0 * (float) pow (n, 4) / 64.0);
    c0 = (15.0 * a * n * n / 16.0) * (1.0 - n +
                                      (3.0 * n * n / 4.0) * (1.0 - n));
    d0 = (35 * a * (float) pow (n, 3) / 48.0) * (1.0 - n +
                                                 11.0 * n * n / 16.0);
    e0 = (315.0 * a * (float) pow (n, 4) / 51.0) * (1.0 - n);

    for (i = 0; i < num_points; i++)
    {
        delta_lon = lon[i] - zone_cm;
        p = delta_lon * PI / 180.0;
        /* convert lat and lon points from decimal degrees to radians */
        lat_rad = lat[i] * PI / 180.0;
//        lon_rad = lon[i]*PI/180.0;

//        rho = a*(1.0-e*e)/(pow((1.0-(e*sin(lat_rad))*(e*sin(lat_rad))),
//              (3.0/2.0)));
        nu = a / pow ((1.0 - (e * sin (lat_rad)) * (e * sin (lat_rad))),
                      (1.0 / 2.0));


        s = a0 * lat_rad - b0 * sin (2 * lat_rad) + c0 * sin (4 * lat_rad) -
            d0 * sin (6 * lat_rad) + e0 * sin (8 * lat_rad);

        /* coefficients for UTM coordinates */
        ki = s * k0;
        kii = nu * sin (lat_rad) * cos (lat_rad) * k0 / 2.0;
        kiii = (pow (nu * sin (lat_rad) * cos (lat_rad), 3) / 24.0) *
            pow ((5 - tan (lat_rad)),
                 2) + 9.0 * eprimesqrd * pow (cos (lat_rad),
                                              2) +
            4.0 * eprimesqrd * eprimesqrd * pow (cos (lat_rad), 4) * k0;
        kiv = nu * cos (lat_rad) * k0;
        kv = pow (cos (lat_rad),
                  3) * (nu / 6.0) * (1 - tan (lat_rad) * tan (lat_rad) +
                                     eprimesqrd * cos (lat_rad) *
                                     cos (lat_rad)) * k0;

        /* calculate UTM coordinates */
        northing = (ki + kii * p * p + kiii * pow (p, 4));
        easting = 500000.0 + (kiv * p + kv * pow (p, 3));

        /* assign to narr_utm array */
        narr_utm[i][0] = easting;
        narr_utm[i][1] = northing;
    }
}


/******************************************************************************
MODULE:  dn_to_radiance

PURPOSE: convert digital counts to radiance for thermal band
         [unit: W m^(-2) sr^(-1) mu^(-1) ]

RETURN: NONE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/30/2014   Song Guo         Original Development
******************************************************************************/
void dn_to_radiance
(
    Input_t *input
)
{
    int pix;
    for (pix = 0; pix < input->size_th.s; pix++)
    {
        if (input->therm_buf[pix] == 0)
        {
            input->therm_buf[pix] = 0;
        }
        else
        {
            input->therm_buf[pix] = (int16_t) (input->meta.gain_th *
                                             (float) input->therm_buf[pix] +
                                             input->meta.bias_th);
        }
        if (input->meta.inst == INST_TM && input->meta.sat == SAT_LANDSAT_5)
        {
            input->therm_buf[pix] = (int16_t) (input->therm_buf[pix] + 0.044);
        }
    }
}


/******************************************************************************
MODULE:  interpolate_to_height

PURPOSE: Inteprolate to height of current pixel

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
10/29/2014   Song Guo
******************************************************************************/
void interpolate_to_height
(
    float *height,        /*I: list of height of one location */
    float *atm1,          /*I: atmospheric parameter colum1 */
    float *atm2,          /*I: atmospheric parameter colum2 */
    float *atm3,          /*I: atmospheric parameter colum3 */
    float interpolate_to, /*I: current landsat pixel height */
    float *at_height      /*O: interpolated height */
)
{
    int i;
    int below = 0;
    int above = 0;
    int close_below;
    int close_above;
    float under_height;
    float under_variables[3];
    float over_height;
    float over_variables[3];
    float m[3];
    float b[3];

    /* Determine points below/above interpolate_to */
    for (i = 0; i < NUM_ELEVATIONS - 1; i++)
    {
        if ((height[i] - interpolate_to) < MINSIGMA)
            below++;
        else
            above++;
    }

    if (below == 0)
        close_below = 0;
    else
        close_below = below;
    under_height = height[close_below];
    under_variables[0] = atm1[close_below];
    under_variables[1] = atm2[close_below];
    under_variables[2] = atm3[close_below];

    if (above == 0)
        close_above = NUM_ELEVATIONS - 1;
    else
        close_below = above;
    over_height = height[close_below];
    over_variables[0] = atm1[close_below];
    over_variables[1] = atm2[close_below];
    over_variables[2] = atm3[close_below];

    for (i = 0; i < 3; i++)
    {
        if (close_above == close_below)
            at_height[i] = under_variables[i];
        else
        {
            m[i] = (over_variables[i] - under_variables[i]) /
                (over_height - under_height);
            b[i] = over_variables[i] - m[i] * over_height;
            at_height[i] = m[i] * interpolate_to + b[i];
        }
    }
}


/******************************************************************************
MODULE:  interpolate_to_location

PURPOSE: Inteprolate to location of current pixel

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
10/29/2014   Song Guo
******************************************************************************/
void interpolate_to_location
(
    float **coordinates,        /*I: current pixel coordinates */
    float **at_height,          /*I: current height atmospheric results */
    float interpolate_easting,  /*I: interpolate to easting */
    float interpolate_northing, /*I: interpolate to northing */
    float *parameters     /*O: interpolated pixel atmospheric parameters */
)
{
    int i, j;
    float h[4];
    float w[4];
    float total = 0.0;

    /* shepard's method */
    for (i = 0; i < 4; i++)
    {
        h[i] = (coordinates[i][0] - interpolate_easting)
                 * (coordinates[i][0] - interpolate_easting)
               + (coordinates[i][1] - interpolate_northing)
                 * (coordinates[i][1] - interpolate_northing);
    }
    qsort (h, 4, sizeof (float), qsort_float_compare_function);

    for (i = 0; i < 4; i++)
    {
        total += 1.0 / h[i];
    }

    for (i = 0; i < 4; i++)
    {
        w[i] = (1.0 / h[i]) / total;
    }

    for (i = 0; i < 3; i++)
    {
        parameters[i] = 0.0;
        for (j = 0; j < 4; j++)
        {
            parameters[i] += w[j] * at_height[i][j];
        }
    }
}


/******************************************************************************
MODULE:  calculate_pixel_atmospheric_parameters

PURPOSE: Generate transmission, upwelled radiance, and downwelled radiance at
         each Landsat pixel

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
10/8/2014   Song Guo         Original Development
******************************************************************************/
int calculate_pixel_atmospheric_parameters
(
    Input_t * input,           /* I: input structure */
    REANALYSIS_POINTS *points, /* I: The coordinate points */
    char *dem_infile,          /* I: address of input DEM filename */
    char *emi_infile,          /* I: address of input Emissivity filename */
    float **modtran_results,   /* I: atmospheric parameter for modtarn run */
    bool verbose               /* I: value to indicate if intermediate
                                     messages be printed */
)
{
    char FUNC_NAME[] = "calculate_pixel_atmospheric_parameters";
    int row;
    int col;
    int line;
    int sample;
    int16_t *dem = NULL;          /* input DEM data (meters) */
    int offset;                 /* offset in the raw binary DEM file to seek to
                                   to begin reading the window in the DEM */
    FILE *dem_fptr = NULL;      /* input scene-based DEM file pointer */
    float **east_grid;
    float **north_grid;
    int first_line;
    float current_easting;
    float current_northing;
    float *distances = NULL;
    int g;
    int n;
    int point;
    int element;
    int closest[6];
//    float easting_near[6];
    float northing_near[6];
    int below[6];
    int min_in_row = 0;
    int min_in_col = 0;
    float **coordinates;
    int indices[4];
    float stay_up, stay_down, stay_right;
    float move_up, move_down, move_right;
    float **at_height;
    int current_index;
    float **current_location;
    float parameters[3];
    char therm_fname[] = "therm_radiance";
    char up_fname[] = "upwell_radiance";
    char down_fname[] = "downwell_radiance";
    char trans_fname[] = "atmos_transmittance";
    FILE *therm_fptr = NULL;
    FILE *trans_fptr = NULL;
    FILE *up_fptr = NULL;
    FILE *down_fptr = NULL;
    float **landsat_results;
    char msg[MAX_STR_LEN];
    char *lst_data_dir = NULL;
    int status;
    int index;

    /* Use local variables for cleaner code */
    int num_rows = points->num_rows;
    int num_cols = points->num_cols;
    int num_points = points->num_points;


    /* Grab the environment path to the LST_DATA_DIR */
    lst_data_dir = getenv ("LST_DATA_DIR");
    if (lst_data_dir == NULL)
    {
        RETURN_ERROR ("LST_DATA_DIR environment variable is not set",
                      FUNC_NAME, FAILURE);
    }

    /* Dynamic allocate the memory */
    east_grid = (float **) allocate_2d_array (num_rows, num_cols,
                                              sizeof (float));
    if (east_grid == NULL)
    {
        RETURN_ERROR ("Allocating east_grid memory", FUNC_NAME, FAILURE);
    }

    north_grid = (float **) allocate_2d_array (num_rows, num_cols,
                                               sizeof (float));
    if (north_grid == NULL)
    {
        RETURN_ERROR ("Allocating north_grid memory", FUNC_NAME, FAILURE);
    }

    index = 0;
    for (row = 0; row < num_rows; row++)
    {
        for (col = 0; col < num_cols; col++)
        {
            east_grid[row][col] = points->utm_easting[index];
            north_grid[row][col] = points->utm_northing[index];

            index++;
        }
    }

    /* Open the DEM for reading raw binary */
    dem_fptr = fopen (dem_infile, "rb");
    if (dem_fptr == NULL)
    {
        RETURN_ERROR ("Error opening the DEM file", FUNC_NAME, FAILURE);
    }

    /* Allocate memory for the DEM */
    dem = (int16_t *) calloc (input->size_th.s, sizeof (int16_t));
    if (dem == NULL)
    {
        RETURN_ERROR ("Error allocating memory for the DEM data",
                      FUNC_NAME, FAILURE);
    }

    /* Open the intermediate binary files for writing
       Note: Needs to be deleted before release */
    therm_fptr = fopen (therm_fname, "wb");
    if (therm_fptr == NULL)
    {
        sprintf (msg, "Opening report file: %s", therm_fname);
        RETURN_ERROR (msg, FUNC_NAME, FAILURE);
    }

    trans_fptr = fopen (trans_fname, "wb");
    if (trans_fptr == NULL)
    {
        sprintf (msg, "Opening report file: %s", trans_fname);
        RETURN_ERROR (msg, FUNC_NAME, FAILURE);
    }

    up_fptr = fopen (up_fname, "wb");
    if (up_fptr == NULL)
    {
        sprintf (msg, "Opening report file: %s", up_fname);
        RETURN_ERROR (msg, FUNC_NAME, FAILURE);
    }

    down_fptr = fopen (down_fname, "wb");
    if (down_fptr == NULL)
    {
        sprintf (msg, "Opening report file: %s", down_fname);
        RETURN_ERROR (msg, FUNC_NAME, FAILURE);
    }

    /* Allocate memory for landsat_results */
    landsat_results = (float **) allocate_2d_array (3, input->size_th.s,
                                                    sizeof (float));
    if (landsat_results == NULL)
    {
        RETURN_ERROR ("Allocating landsat_results memory", FUNC_NAME, FAILURE);
    }

    /* Allocate memory for coordinates */
    coordinates = (float **) allocate_2d_array (4, 2, sizeof (float));
    if (coordinates == NULL)
    {
        RETURN_ERROR ("Allocating coordinates memory", FUNC_NAME, FAILURE);
    }

    /* Allocate memory for at_height */
    at_height = (float **) allocate_2d_array (4, 3, sizeof (float));
    if (at_height == NULL)
    {
        RETURN_ERROR ("Allocating at_height memory", FUNC_NAME, FAILURE);
    }

    /* Allocate memory for current_location */
    current_location =
        (float **) allocate_2d_array (num_points * NUM_ELEVATIONS,
                                      LST_NUM_ELEMENTS, sizeof (float));
    if (current_location == NULL)
    {
        RETURN_ERROR ("Allocating current_location memory", FUNC_NAME,
                      FAILURE);
    }

    if (verbose)
    {
        LOG_MESSAGE ("Iterate through all lines in landsat scene\n",
                     FUNC_NAME);
    }

    /* Loop through each line in the image */
    for (line = 0; line < input->size_th.l; line++)
    {
        /* Print status on every 1000 lines */
        if (!(line % 250))
        {
            if (verbose)
            {
                printf ("Processing line %d\r", line);
                fflush (stdout);
            }
        }

        /* Read the input thermal band -- data is read into input->therm_buf */
        if (!GetInputThermLine (input, line))
        {
            sprintf (msg, "Reading input thermal data for line %d", line);
            RETURN_ERROR (msg, FUNC_NAME, FAILURE);
        }
        dn_to_radiance (input);

        /* Can also read in one line of DEM data here */
        /* Start reading DEM from the start_line */
        offset = sizeof (int16_t) * line * input->size_th.s;
        fseek (dem_fptr, offset, SEEK_SET);
        if (fread (dem, sizeof (int16_t), input->size_th.s, dem_fptr)
            != input->size_th.s)
        {
            sprintf (msg, "Error reading values from the DEM file "
                             "starting at line %d.", line);
            RETURN_ERROR (msg, FUNC_NAME, FAILURE);
        }

        /* set first_line be 1 */
        first_line = 1;
        for (sample = 0; sample < input->size_th.s; sample++)
        {
// TODO TODO TODO - Change to use fill value from xml, not 0
// TODO TODO TODO - Change to use fill value from xml, not 0
// TODO TODO TODO - Change to use fill value from xml, not 0
// TODO TODO TODO - Change to use fill value from xml, not 0
// TODO TODO TODO - Change to use fill value from xml, not 0
// TODO TODO TODO - Change to use fill value from xml, not 0
// TODO TODO TODO - Change to use fill value from xml, not 0
// TODO TODO TODO - Change to use fill value from xml, not 0
// TODO TODO TODO - Change to use fill value from xml, not 0
// TODO TODO TODO - Change to use fill value from xml, not 0
// TODO TODO TODO - Change to use fill value from xml, not 0
// TODO TODO TODO - Change to use fill value from xml, not 0
// TODO TODO TODO - Change to use fill value from xml, not 0
            if (input->therm_buf != 0)
            {
                /* determine UTM coordinates of current pixel */
                current_easting = input->meta.ul_map_corner.x
                    + line * input->meta.pixel_size[0];
                current_northing = input->meta.ul_map_corner.y
                    + line * input->meta.pixel_size[1];

                if (first_line == 1)
                {
                    /* compute distance between current pixel and each narr
                       points in UTM coordinates

                       Note: consider only calculating points within a small
                       range nearby */
                    distances = malloc (num_points * sizeof (float));
                    if (distances == NULL)
                    {
                        ERROR_MESSAGE ("Allocating distances memory",
                                       FUNC_NAME);
                    }
                    for (point = 0; point < num_points; point++)
                    {
                        distance_in_utm (points->utm_easting[point],
                                         points->utm_northing[point],
                                         current_easting,
                                         current_northing,
                                         &distances[point]);
                    }

                    /* find indices of 6 closet points */
LOG_MESSAGE ("Before qsort", FUNC_NAME);
                    qsort (distances, num_points, sizeof (float),
                           qsort_float_compare_function);
LOG_MESSAGE ("After qsort", FUNC_NAME);

                    n = 0;
                    /* find indices of 6 closest points */
                    for (g = 0; g < 6; g++)
                    {
                        for (point = 0; point < num_points; point++)
                        {
                            {
                                closest[g] = point;
//                                easting_near[g] = points->utm_easting[g];
                                northing_near[g] = points->utm_northing[g];
                                if ((northing_near[g] - current_northing) <=
                                    MINSIGMA)
                                {
                                    below[n] = point;
                                    n++;
                                }
                            }
                        }
                    }
LOG_MESSAGE ("HERE 1", FUNC_NAME);

                    min_in_row = points->row[closest[below[0]]];
                    min_in_col = points->col[closest[below[0]]];
                    for (g = 1; g < n; g++)
                    {
                        min_in_row =
                            min (min_in_row, points->row[closest[below[g]]]);
                        min_in_col =
                            min (min_in_col, points->col[closest[below[g]]]);
                    }
LOG_MESSAGE ("HERE 2", FUNC_NAME);

                    /* extract UTM coordinates of four points to be
                       interpolated and build array */
                    coordinates[0][0] = east_grid[min_in_row][min_in_col];
                    coordinates[0][1] = north_grid[min_in_row][min_in_col];
                    coordinates[1][0] = east_grid[min_in_row][min_in_col + 1];
                    coordinates[1][1] = north_grid[min_in_row][min_in_col + 1];
                    coordinates[2][0] = east_grid[min_in_row + 1][min_in_col];
                    coordinates[2][1] = north_grid[min_in_row + 1][min_in_col];
                    coordinates[3][0] =
                        east_grid[min_in_row + 1][min_in_col + 1];
                    coordinates[3][1] =
                        north_grid[min_in_row + 1][min_in_col + 1];

                    /* determine index of four points in order to pull from
                       MODTRAN results */
                    indices[0] = min_in_row * num_cols + min_in_col;
                    indices[1] = (min_in_row + 1) * num_cols + min_in_col;
                    indices[2] = min_in_row * num_cols + min_in_col + 1;
                    indices[3] = (min_in_row + 1) * num_cols + min_in_col + 1;

                    /* set firstInLine variable to false */
                    first_line = 0;

                    free (distances);
                    distances = NULL;
LOG_MESSAGE ("HERE XXX", FUNC_NAME);
                }
                else
                {
LOG_MESSAGE ("****** NOT ****** First Line", FUNC_NAME);
                    /* given indices of previous pixel, there are six possible
                       quads to move into check 6 distances to determine new
                       upperleft corner */
                    distance_in_utm (east_grid[min_in_row][min_in_col],
                                     north_grid[min_in_row][min_in_col],
                                     current_easting, current_northing,
                                     &stay_right);

                    if ((min_in_col + 2) < num_cols)
                        distance_in_utm (east_grid[min_in_row][min_in_col + 2],
                                         north_grid[min_in_row][min_in_col + 2],
                                         current_easting, current_northing,
                                         &move_right);
                    else
                        stay_right = (float) SHRT_MAX;

                    distance_in_utm (east_grid[min_in_row + 1][min_in_col + 1],
                                     north_grid[min_in_row + 1][min_in_col + 1],
                                     current_easting, current_northing,
                                     &stay_up);

                    distance_in_utm (east_grid[min_in_row - 1][min_in_col + 1],
                                     north_grid[min_in_row - 1][min_in_col + 1],
                                     current_easting, current_northing,
                                     &move_up);

                    distance_in_utm (east_grid[min_in_row][min_in_col + 1],
                                     north_grid[min_in_row][min_in_col + 1],
                                     current_easting, current_northing,
                                     &stay_down);

                    if ((min_in_row + 2) < num_rows)
                    {
                        distance_in_utm (
                            east_grid[min_in_row + 2][min_in_col + 1],
                            north_grid[min_in_row + 2][min_in_col + 1],
                            current_easting, current_northing, &move_down);
                    }
                    else
                        move_down = (float) SHRT_MAX;

                    if ((move_right - stay_right) < MINSIGMA)
                        min_in_col++;
                    if ((move_up - stay_up) < MINSIGMA)
                        min_in_row--;
                    if ((move_down - stay_down) < MINSIGMA)
                        min_in_row++;

                    /* extract UTM coordinates of four points to be interpolated
                       and build array */
                    coordinates[0][0] = east_grid[min_in_row][min_in_col];
                    coordinates[0][1] = north_grid[min_in_row][min_in_col];
                    coordinates[1][0] = east_grid[min_in_row + 1][min_in_col];
                    coordinates[1][1] = north_grid[min_in_row + 1][min_in_col];
                    coordinates[2][0] = east_grid[min_in_row][min_in_col + 1];
                    coordinates[2][1] = north_grid[min_in_row][min_in_col + 1];
                    coordinates[3][0] =
                        east_grid[min_in_row + 1][min_in_col + 1];
                    coordinates[3][1] =
                        north_grid[min_in_row + 1][min_in_col + 1];

                    /* determine index of four points in order to pull from
                       MODTRAN results */
                    indices[0] = min_in_row * num_cols + min_in_col;
                    indices[1] = (min_in_row + 1) * num_cols + min_in_col;
                    indices[2] = min_in_row * num_cols + min_in_col + 1;
                    indices[3] = (min_in_row + 1) * num_cols + min_in_col + 1;
                }

                /* convert height from m to km */
                dem[sample] = (float) dem[sample] / 1000.0;

                /* interpolate three parameters to that height at each of the
                   four closest points */
                for (g = 0; g < 4; g++)
                {
                    current_index = indices[g] * NUM_ELEVATIONS;
                    /* extract atmospheric parameters for all heights at current
                       location */
                    for (element = 0; element < LST_NUM_ELEMENTS; element++)
                    {
                        for (index = current_index;
                             index < current_index + NUM_ELEVATIONS - 1;
                             index++)
                        {
                            current_location[index - current_index][element] =
                                modtran_results[index][element];
                        }
                    }

                    /* interpolate three atmospheric parameters to current
                       height */
                    interpolate_to_height (current_location[LST_HEIGHT],
                        current_location[LST_TRANSMISSION],
                        current_location[LST_UPWELLED_RADIANCE],
                        current_location[LST_DOWNWELLED_RADIANCE],
                        dem[sample], at_height[g]);
                }

                /* interpolate parameters at appropriate height to location of
                   current pixel */
                interpolate_to_location (coordinates, at_height,
                                         current_easting, current_northing,
                                         parameters);

                /* convert radiances to W*m^(-2)*sr(-1) */
                landsat_results[0][sample] = parameters[0];
                landsat_results[1][sample] = 10000.0 * parameters[1];
                landsat_results[2][sample] = 10000.0 * parameters[2];
            } /* END - if (input->therm_buf != 0) */
        } /* END - for sample */

        /* Write out the temporary binary output files
           Note: needs to be deleted before release */
        status = fwrite (&input->therm_buf[0], sizeof (int16_t),
                         input->size_th.s, therm_fptr);
        if (status != input->size_th.s)
        {
            sprintf (msg, "Writing to %s", therm_fname);
            ERROR_MESSAGE (msg, FUNC_NAME);
        }

        status = fwrite (&landsat_results[0][0], sizeof (float),
                         input->size_th.s, trans_fptr);
        if (status != input->size_th.s)
        {
            sprintf (msg, "Writing to %s", trans_fname);
            ERROR_MESSAGE (msg, FUNC_NAME);
        }

        status = fwrite (&landsat_results[1][0], sizeof (float),
                         input->size_th.s, up_fptr);
        if (status != input->size_th.s)
        {
            sprintf (msg, "Writing to %s", up_fname);
            ERROR_MESSAGE (msg, FUNC_NAME);
        }

        status = fwrite (&landsat_results[2][0], sizeof (float),
                         input->size_th.s, down_fptr);
        if (status != input->size_th.s)
        {
            sprintf (msg, "Writing to %s", down_fname);
            ERROR_MESSAGE (msg, FUNC_NAME);
        }
    } /* END - for line */

    /* Free allocated memory */
    status = free_2d_array ((void **) current_location);
    if (status != SUCCESS)
    {
        ERROR_MESSAGE ("Freeing memory: at_height\n", FUNC_NAME);
    }
    status = free_2d_array ((void **) at_height);
    if (status != SUCCESS)
    {
        ERROR_MESSAGE ("Freeing memory: at_height\n", FUNC_NAME);
    }
    status = free_2d_array ((void **) coordinates);
    if (status != SUCCESS)
    {
        ERROR_MESSAGE ("Freeing memory: coordinates\n", FUNC_NAME);
    }
    status = free_2d_array ((void **) east_grid);
    if (status != SUCCESS)
    {
        ERROR_MESSAGE ("Freeing memory: east_grid\n", FUNC_NAME);
    }
    status = free_2d_array ((void **) north_grid);
    if (status != SUCCESS)
    {
        ERROR_MESSAGE ("Freeing memory: north_grid\n", FUNC_NAME);
    }

    status = free_2d_array ((void **) landsat_results);
    if (status != SUCCESS)
    {
        ERROR_MESSAGE ("Freeing memory: landsat_results\n", FUNC_NAME);
    }

    /* Close the intermediate binary files
       Note: needs to be deleted before release */
    status = fclose (therm_fptr);
    if (status)
    {
        sprintf (msg, "Closing file %s", therm_fname);
        ERROR_MESSAGE (msg, FUNC_NAME);
    }

    status = fclose (trans_fptr);
    if (status)
    {
        sprintf (msg, "Closing file %s", trans_fname);
        ERROR_MESSAGE (msg, FUNC_NAME);
    }

    status = fclose (up_fptr);
    if (status)
    {
        sprintf (msg, "Closing file %s", up_fname);
        ERROR_MESSAGE (msg, FUNC_NAME);
    }

    status = fclose (down_fptr);
    if (status)
    {
        sprintf (msg, "Closing file %s", down_fname);
        ERROR_MESSAGE (msg, FUNC_NAME);
    }

    return SUCCESS;
}
