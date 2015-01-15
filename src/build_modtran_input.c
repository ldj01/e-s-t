
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>


#include "const.h"
#include "2d_array.h"
#include "utilities.h"
#include "input.h"
#include "scene_based_lst.h"


#define STANRDARD_GRAVITY_IN_M_PER_SEC_SQRD 9.80665
#define POLAR_RADIUS_IN_KM 6356.752
#define EQUATORIAL_RADIUS_IN_KM 6378.137


/*****************************************************************************
MODULE:  convert_geopotential_geometric

PURPOSE: Convert array of geopotential heights to array of geometric heights
         given latitude

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/22/2014   Song Guo         Original Development
*****************************************************************************/
int convert_geopotential_geometric
(
    int num_points,        /* I: number of points */
    float *lat,            /* I: latitude in degrees */
    float **geo_potential, /* I: geo_potential height */
    float **geo_metric     /* O: geo_metric height */
)
{
    char FUNC_NAME[] = "convert_geopotential_geometric";
    int i, j;
    double *radlat;        /* radians of the latitude */
    double *radius;        /* radius in meters at the latitude */
    double *gravity_ratio; /* gravity ratio at the latitude */

    double sin_lat;
    double cos_lat;
    double cos_2lat;

    double inv_r_min_sqrd;
    double inv_r_max_sqrd;
    double inv_std_gravity;

    int array_bytes = num_points * sizeof (double);

    /* Allocate memeory */
    radlat = (double *) malloc (array_bytes);
    if (radlat == NULL)
    {
        RETURN_ERROR ("Allocating radlat memory", FUNC_NAME, FAILURE);
    }

    radius = (double *) malloc (array_bytes);
    if (radius == NULL)
    {
        RETURN_ERROR ("Allocating radius memory", FUNC_NAME, FAILURE);
    }

    gravity_ratio = (double *) malloc (array_bytes);
    if (gravity_ratio == NULL)
    {
        RETURN_ERROR ("Allocating gravity_ratio memory", FUNC_NAME, FAILURE);
    }

    inv_r_min_sqrd = 1.0 / (POLAR_RADIUS_IN_KM * POLAR_RADIUS_IN_KM);
    inv_r_max_sqrd = 1.0 / (EQUATORIAL_RADIUS_IN_KM * EQUATORIAL_RADIUS_IN_KM);
    inv_std_gravity = 1.0 / STANRDARD_GRAVITY_IN_M_PER_SEC_SQRD;

    /* Calculate values for each point
       NOTE - These are constant for this scene, yet calculated twice since
              this routine is called twice.  Should not be a time hog though
              since the number of points is low and will be lower with
              MERRA data if we use MERRA later. */
    for (i = 0; i < num_points; i++)
    {
        /* Convert the latitude to radians */
        radlat[i] = lat[i] * RAD;

        /* Only calculate these once for the latitude */
        sin_lat = sin (radlat[i]);
        cos_lat = cos (radlat[i]);
        cos_2lat = cos (2.0 * radlat[i]);

        /* Determine the radius at the latitude in meters */
        /* TODO TODO TODO - if the r_min and r_max were meters to start with,
                            then we would not have to multiply by 1000.0 here
                            to get meters.  Or at least that is what I
                            currently assume is happening. */
        /* TODO TODO TODO - if the r_min and r_max were meters to start with,
                            then we would not have to multiply by 1000.0 here
                            to get meters.  Or at least that is what I
                            currently assume is happening. */
        /* TODO TODO TODO - if the r_min and r_max were meters to start with,
                            then we would not have to multiply by 1000.0 here
                            to get meters.  Or at least that is what I
                            currently assume is happening. */
        radius[i] = 1000.0
                    * sqrt (1.0 / (((cos_lat * cos_lat) * inv_r_max_sqrd)
                                   + ((sin_lat * sin_lat) * inv_r_min_sqrd)));

        /* Determine the gravity ratio at the latitude */
        gravity_ratio[i] = (9.80616
                            * (1.0
                               - (0.002637 * cos_2lat)
                               + (0.0000059 * (cos_2lat * cos_2lat))))
                           * inv_std_gravity;
    }

    /* Calculate the geompetric height for each point in each layer */
    for (i = 0; i < P_LAYER; i++)
    {
        for (j = 0; j < num_points; j++)
        {
            /* The values should all be double until converted to float here */
            geo_metric[i][j] = (geo_potential[i][j] * radius[j])
                               / (1000.0 * (gravity_ratio[j] * radius[j]
                                            - geo_potential[i][j]));
        }
    }

    /* free memory */
    free (radlat);
    free (radius);
    free (gravity_ratio);

    return SUCCESS;
}


/*****************************************************************************
MODULE:  convert_sh_rh

PURPOSE: Given array of specific humidities, temperature, and pressure,
         generate array of relative humidities

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/22/2014   Song Guo         Original Development
*****************************************************************************/
int convert_sh_rh
(
    int num_points,   /* I: number of points */
    float *lat,       /* I: latitude in degrees */
    float **spec_hum,
    float **temp_k,
    float **pressure,
    float **rh        /* O: relative humidity */
)
{
    char FUNC_NAME[] = "convert_sh_rh";
    int i, j;
    float mh20 = 18.01534;
    float mdry = 28.9644;

    float a0w = 6.107799961;
    float a1w = 4.436518521e-1;
    float a2w = 1.428945805e-2;
    float a3w = 2.650648471e-4;
    float a4w = 3.0312403963e-6;
    float a5w = 2.034080948e-8;
    float a6w = 6.136820929e-11;

    float **temp_c;
    float **ewater;
// TODO TODO TODO - RDD - This was calculated but not used??????????
//    float **e2;
    float **goff;
    float **ph20;

    /* Allocate memory */
    temp_c = (float **) allocate_2d_array (P_LAYER, num_points,
                                           sizeof (float));
    if (temp_c == NULL)
    {
        RETURN_ERROR ("Allocating temp_c memory", FUNC_NAME, FAILURE);
    }

    ewater = (float **) allocate_2d_array (P_LAYER, num_points,
                                           sizeof (float));
    if (ewater == NULL)
    {
        RETURN_ERROR ("Allocating ewater memory", FUNC_NAME, FAILURE);
    }

// TODO TODO TODO - RDD - This was calculated but not used??????????
//    e2 = (float **) allocate_2d_array (P_LAYER, num_points, sizeof (float));
//    if (e2 == NULL)
//    {
//        RETURN_ERROR ("Allocating e2 memory", FUNC_NAME, FAILURE);
//    }

    goff = (float **) allocate_2d_array (P_LAYER, num_points, sizeof (float));
    if (goff == NULL)
    {
        RETURN_ERROR ("Allocating goff memory", FUNC_NAME, FAILURE);
    }

    ph20 = (float **) allocate_2d_array (P_LAYER, num_points, sizeof (float));
    if (ph20 == NULL)
    {
        RETURN_ERROR ("Allocating  memory", FUNC_NAME, FAILURE);
    }

    for (i = 0; i < P_LAYER; i++)
    {
        for (j = 0; j < num_points; j++)
        {
            /* Convert temperature to C */
            temp_c[i][j] = temp_k[i][j] - 273.15;

            /* calculate vapor pressure at given temperature */
            ewater[i][j] = a0w + temp_c[i][j]
                           * (a1w + temp_c[i][j]
                              * (a2w + temp_c[i][j]
                                 * (a3w + temp_c[i][j]
                                    * (a4w + temp_c[i][j]
                                       * (a5w + temp_c[i][j]
                                          * (a6w * temp_c[i][j]))))));/* hpa */

// TODO TODO TODO - RDD - This was calculated but not used??????????
//            e2[i][j] = exp (-0.58002206e4 / temp_k[i][j]
//                            + 0.13914993
//                            - 0.48640239e-1 * temp_k[i][j]
//                            + 0.41764768e-4 * pow (temp_k[i][j], 2.0)
//                            - 0.14452093e-7 * pow (temp_k[i][j], 3.0)
//                            + 0.65459673 * log (temp_k[i][j])); /* Pa */

            goff[i][j] = -7.90298 * (373.16 / temp_k[i][j] - 1.0)
                         + 5.02808 * log10 (373.16 / temp_k[i][j])
                         - 1.3816e-7 * pow (10.0, (11.344
                                                   * (1.0
                                                      - (temp_k[i][j]
                                                         / 373.16)))
                                                  - 1.0)
                         + 8.1328e-3 * pow (10.0, (-3.49149
                                                   * (373.16 / temp_k[i][j]
                                                      - 1.0))
                                                  - 1.0)
                         + log10 (1013.246); /* hPa */

            ph20[i][j] = (spec_hum[i][j] * pressure[i][j] * mdry)
                         / (mh20
                            - spec_hum[i][j] * mh20
                            + spec_hum[i][j] * mdry);

            rh[i][j] = (ph20[i][j] / pow (10.0, goff[i][j])) * 100.0;
        }
    }

    /* Free allocated memory */
    if (free_2d_array ((void **) temp_c) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: temp_c\n", FUNC_NAME, FAILURE);
    }

    if (free_2d_array ((void **) ewater) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: ewater\n", FUNC_NAME, FAILURE);
    }

// TODO TODO TODO - RDD - This was calculated but not used??????????
//    if (free_2d_array ((void **) e2) != SUCCESS)
//    {
//        RETURN_ERROR ("Freeing memory: e2\n", FUNC_NAME, FAILURE);
//    }

    if (free_2d_array ((void **) goff) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: goff\n", FUNC_NAME, FAILURE);
    }

    if (free_2d_array ((void **) ph20) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: ph20\n", FUNC_NAME, FAILURE);
    }

    return SUCCESS;
}


/******************************************************************************
MODULE:  read_narr_coordinates

PURPOSE: Creates directories and writes tape5 file, caseList, and commandList

RETURN: SUCCESS
        FAILURE
******************************************************************************/
int read_narr_coordinates
(
    char *lst_data_dir,
    int **grid_i,
    int **grid_j,
    float **lat,
    float **lon
)
{
    char FUNC_NAME[] = "read_narr_coordinates";
    int i, j;
    int count;
    char coord_file[PATH_MAX];
    FILE *fd = NULL;

    /* Setup the string to be used to open the coordinates file */
    count = snprintf (coord_file, sizeof (coord_file),
                      "%s/%s", lst_data_dir, "narr_coordinates.txt");
    if (count < 0 || count >= sizeof (coord_file))
    {
        RETURN_ERROR ("Failed initializing coord_file variable for"
                      " narr_coordinates.txt", FUNC_NAME, FAILURE);
    }

    /* Open the NARR coordinates file */
    fd = fopen (coord_file, "r");
    if (fd == NULL)
    {
        RETURN_ERROR ("Can't open narr_coordinates.txt file", FUNC_NAME,
                      FAILURE);
    }

    /* Read each value from the file */
    for (j = 0; j < NARR_COL; j++)
    {
        for (i = 0; i < NARR_ROW; i++)
        {
            /* File Format:
               Grid_I Grid_J Latitude Longitude
             */
            if (fscanf (fd, "%d %d %f %f",
                        &grid_i[i][j], &grid_j[i][j],
                        &lat[i][j], &lon[i][j]) == EOF)
            {
                RETURN_ERROR ("End of file (EOF) is met before"
                              " NARR_ROW * NARR_COL lines",
                              FUNC_NAME, FAILURE);
            }

            /* TODO - Should think about fixing the input file, so that this
                      confusing conversion is not needed.
                      When you read  the file data, it is as if you are
                      reading the values from the lower left to the upper
                      right as applied to the earth.  And the values being
                      read in start with an origin somewhere around the lower
                      left, hence the need for the following conversion. */
            if ((lon[i][j] - 180.0) > MINSIGMA)
                lon[i][j] = 360.0 - lon[i][j];
            else
                lon[i][j] = -lon[i][j];
        }
    }

    if (fclose (fd) != SUCCESS)
    {
        RETURN_ERROR ("Closing file: narr_coordinates.txt\n",
                      FUNC_NAME, FAILURE);
    }

    return SUCCESS;
}


/******************************************************************************
MODULE:  read_std_mid_lat_summer_atmos

PURPOSE: Creates directories and writes tape5 file, caseList, and commandList

RETURN: SUCCESS
        FAILURE
******************************************************************************/
int read_std_mid_lat_summer_atmos
(
    char *lst_data_dir,
    float *stan_height,
    float *stan_pre,
    float *stan_temp,
    float *stan_rh
)
{
    char FUNC_NAME[] = "read_std_mid_lat_summer_atmos";
    int i;
    int count;
    char atmos_file[PATH_MAX];
    FILE *fd = NULL;

    /* read in file containing standard mid lat summer atmosphere information
       to be used for upper layers */
    count = snprintf (atmos_file, sizeof (atmos_file),
                      "%s/%s", lst_data_dir, "std_mid_lat_summer_atmos.txt");
    if (count < 0 || count >= sizeof (atmos_file))
    {
        RETURN_ERROR ("Failed initializing atmos_file variable for"
                      " std_mid_lat_summer_atmos.txt", FUNC_NAME, FAILURE);
    }

    fd = fopen (atmos_file, "r");
    if (fd == NULL)
    {
        RETURN_ERROR ("Opening file: std_mid_lat_summer_atmos.txt\n",
                      FUNC_NAME, FAILURE);
    }

    for (i = 0; i < STAN_LAYER; i++)
    {
        if (fscanf (fd, "%f %f %f %f", &stan_height[i], &stan_pre[i],
                    &stan_temp[i], &stan_rh[i]) == EOF)
        {
            RETURN_ERROR ("End of file (EOF) is met before STAN_LAYER lines",
                          FUNC_NAME, FAILURE);
        }
    }

    if (fclose (fd) != SUCCESS)
    {
        RETURN_ERROR ("Closing file: std_mid_lat_summer_atmos.txt\n",
                      FUNC_NAME, FAILURE);
    }

    return SUCCESS;
}


/******************************************************************************
MODULE:  read_narr_parameter_values

PURPOSE: Creates directories and writes tape5 file, caseList, and commandList

RETURN: SUCCESS
        FAILURE
******************************************************************************/
int read_narr_parameter_values
(
    int *layers,
    char *parameter,
    float **output_2d_array
)
{
    char FUNC_NAME[] = "read_narr_parameter_values";
    char msg_str[MAX_STR_LEN];
    char parm_filename[PATH_MAX];
    int i, j;
    int count;
    int file_rows, file_cols;
    FILE *fd = NULL;

    /* Read each layers parameter file into memory */
    for (i = 0; i < P_LAYER; i++)
    {
        /* Build the full path to the parameter file */
        count = snprintf (parm_filename, sizeof (parm_filename),
                           "%s/%d%s", parameter, layers[i], ".txt");
        if (count < 0 || count >= sizeof (parm_filename))
        {
            RETURN_ERROR ("Failed initializing parm_filename variable for"
                          " NARR parameters", FUNC_NAME, FAILURE);
        }

        snprintf(msg_str, sizeof (msg_str), "Reading [%s]", parm_filename);
        LOG_MESSAGE (msg_str, FUNC_NAME);

        /* Open the parameter file for reading */
        fd = fopen (parm_filename, "r");
        if (fd == NULL)
        {
            RETURN_ERROR ("Can't HGT_1 txt file", FUNC_NAME, FAILURE);
        }

        /* Read the number of rows and columns and verify them */
        fscanf (fd, "%d %d", &file_rows, &file_cols);
        if (file_rows != NARR_ROW || file_cols != NARR_COL)
        {
            RETURN_ERROR ("Parameter file contains invalid number of rows and"
                          " columns", FUNC_NAME, FAILURE);
        }

        /* Read the values into memory */
        for (j = 0; j < NARR_ROW * NARR_COL; j++)
        {
            if (fscanf (fd, "%f", &output_2d_array[i][j]) == EOF)
            {
                RETURN_ERROR ("End of file (EOF) is met before "
                              "NARR_ROW * NARR_COL lines",
                               FUNC_NAME, FAILURE);
            }
        }

        /* Close the current parameter file */
        if (fclose (fd) != SUCCESS)
        {
            snprintf(msg_str, sizeof (msg_str),
                     "Closing file: %s/%d.txt\n", parameter, layers[i]);
            RETURN_ERROR (msg_str, FUNC_NAME, FAILURE);
        }
    }

    return SUCCESS;
}


/******************************************************************************
MODULE:  build_modtran_input

PURPOSE: Creates directories and writes tape5 file, caseList, and commandList

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/21/2014   Song Guo         Original Development
******************************************************************************/
int build_modtran_input
(
    Input_t *input,  /*I: input structure */
    int *num_pts,    /*O: number of NARR points */
    bool verbose     /*I: value to indicate if intermediate messages should
                          be printed */
)
{
    char FUNC_NAME[] = "build_modtran_input";
    int **eye;
    int **jay;
    float **lat;
    float **lon;
    float **hgt1;
    float **spfh1;
    float **tmp1;
    float **hgt2;
    float **spfh2;
    float **tmp2;
    float *narr_lat;
    float *narr_lon;
    float **narr_hgt1;
    float **narr_spfh1;
    float **narr_tmp1;
    float **narr_hgt2;
    float **narr_spfh2;
    float **narr_tmp2;
    float **narr_height;
    float **narr_height1;
    float **narr_height2;
    float **narr_rh;
    float **narr_rh1;
    float **narr_rh2;
    float **narr_tmp;
    int i, j, k;
    int layers[P_LAYER] = { 1000, 975, 950, 925, 900,
                            875, 850, 825, 800, 775,
                            750, 725, 700, 650, 600,
                            550, 500, 450, 400, 350,
                            300, 275, 250, 225, 200,
                            175, 150, 125, 100 };
    int num_points;
    float narr_ul_lat;
    float narr_ul_lon;
    float narr_lr_lat;
    float narr_lr_lon;
    int max_eye;
    int min_eye;
    int max_jay;
    int min_jay;
    int num_eyes;
    int num_jays;
    float **pressure;
    int rem1;
    int rem2;
    float hour1;
    float hour2;
    float inv_hour_diff;
    float time;
    float time_diff;
    FILE *fd;
    float *stan_height;
    float *stan_pre;
    float *stan_temp;
    float *stan_rh;
    float *temp_height;
    float *temp_pressure;
    float *temp_temp;
    float *temp_rh;
    float gndalt[NUM_ELEVATIONS] = { 0.0, 0.6, 1.1, 1.6, 2.1, 2.6,
        3.1, 3.6, 4.05
    };
    float inv_narr_height_diff;
    int num_modtran_runs;
    char command[MAX_STR_LEN];
    char current_gdalt[MAX_STR_LEN];
    char current_temp[MAX_STR_LEN];
    char current_alb[MAX_STR_LEN];
    char current_point[MAX_STR_LEN];
    char temp_out[MAX_STR_LEN];
    char curr_path[MAX_STR_LEN];
    int index_below = 0;
    int index_above = NUM_ELEVATIONS;
    float new_height;
    float new_pressure;
    float new_temp;
    float new_rh;
    int index, index2;
    int *counter;
    float tmp[3] = { 273.0, 310.0, 0.0 };
    float alb[3] = { 0.0, 0.0, 0.1 };
    char *lst_data_dir = NULL;
    char *modtran_path = NULL;
    char *modtran_data_dir = NULL;
    int case_counter;
    char **case_list;
    char **command_list;
    char lat_str[7]; /* 6 plus the string termination character */
    char lon_str[7]; /* 6 plus the string termination character */
    char msg_str[MAX_STR_LEN];

    /* Dynamic allocate the 2d memory for the coordinates */
    eye = (int **) allocate_2d_array (NARR_ROW, NARR_COL, sizeof (int));
    if (eye == NULL)
    {
        RETURN_ERROR ("Allocating eye memory", FUNC_NAME, FAILURE);
    }

    jay = (int **) allocate_2d_array (NARR_ROW, NARR_COL, sizeof (int));
    if (jay == NULL)
    {
        RETURN_ERROR ("Allocating jay memory", FUNC_NAME, FAILURE);
    }

    lat = (float **) allocate_2d_array (NARR_ROW, NARR_COL, sizeof (float));
    if (lat == NULL)
    {
        RETURN_ERROR ("Allocating lat memory", FUNC_NAME, FAILURE);
    }

    lon = (float **) allocate_2d_array (NARR_ROW, NARR_COL, sizeof (float));
    if (lon == NULL)
    {
        RETURN_ERROR ("Allocating lon memory", FUNC_NAME, FAILURE);
    }

    /* Grab the environment path to the LST_DATA_DIR */
    lst_data_dir = getenv ("LST_DATA_DIR");
    if (lst_data_dir == NULL)
    {
        RETURN_ERROR ("LST_DATA_DIR environment variable is not set",
                      FUNC_NAME, FAILURE);
    }

    /* Read the coordinates into memory */
    if (read_narr_coordinates (lst_data_dir, eye, jay, lat, lon) != SUCCESS)
    {
        RETURN_ERROR ("Failed loading HGT_1 parameters", FUNC_NAME, FAILURE);
    }

    /* expand range to include NARR points outside image for edge pixels */
    narr_ul_lat = input->meta.ul_geo_corner.lat + 0.5;
    narr_ul_lon = input->meta.ul_geo_corner.lon - 0.5;
    narr_lr_lat = input->meta.lr_geo_corner.lat - 0.5;
    narr_lr_lon = input->meta.lr_geo_corner.lon + 0.5;

    /* determine what points in the NARR dataset fall within the Landsat
       image using logical operators lessThanLat and greaterThanLat are values
       where the NARR values are less than or greater than the edges of the
       Landsat corners values respectively pixels that are true in both fall
       within the Landsat scene the same thing is done with longitude values */
    max_eye = 0;
    min_eye = 1000;
    max_jay = 0;
    min_jay = 1000;
    for (i = 0; i < NARR_ROW - 1; i++)
    {
        for (j = 0; j < NARR_COL - 1; j++)
        {
            if ((lat[i][j] - narr_ul_lat) < MINSIGMA
                && (lat[i][j] - narr_lr_lat) > MINSIGMA
                && (lon[i][j] - narr_lr_lon) < MINSIGMA
                && (lon[i][j] - narr_ul_lon) > MINSIGMA)
            {
                max_eye = max (max_eye, eye[i][j]);
                min_eye = min (min_eye, eye[i][j]);
                max_jay = max (max_jay, jay[i][j]);
                min_jay = min (min_jay, jay[i][j]);
            }
        }
    }
    max_eye--;
    min_eye--;
    max_jay--;
    min_jay--;
    num_eyes = max_eye - min_eye + 1;
    num_jays = max_jay - min_jay + 1;
    num_points = num_eyes * num_jays;

    /* Free eye and jay memory since they are not used anymore */
    if (free_2d_array ((void **) eye) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: eye\n", FUNC_NAME, FAILURE);
    }
    eye = NULL;

    if (free_2d_array ((void **) jay) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: jay\n", FUNC_NAME, FAILURE);
    }
    jay = NULL;

    /* ==================================================================== */

    /* Dynamic allocate the 2d memory for time before Landsat Acq. */
    hgt1 = (float **) allocate_2d_array (P_LAYER, NARR_ROW * NARR_COL,
                                         sizeof (float));
    if (hgt1 == NULL)
    {
        RETURN_ERROR ("Allocating HGT_1 memory", FUNC_NAME, FAILURE);
    }

    spfh1 = (float **) allocate_2d_array (P_LAYER, NARR_ROW * NARR_COL,
                                          sizeof (float));
    if (spfh1 == NULL)
    {
        RETURN_ERROR ("Allocating SPFH_1 memory", FUNC_NAME, FAILURE);
    }

    tmp1 = (float **) allocate_2d_array (P_LAYER, NARR_ROW * NARR_COL,
                                         sizeof (float));
    if (tmp1 == NULL)
    {
        RETURN_ERROR ("Allocating TMP_1 memory", FUNC_NAME, FAILURE);
    }

    /* Dynamic allocate the 2d memory for time after Landsat Acq. */
    hgt2 = (float **) allocate_2d_array (P_LAYER, NARR_ROW * NARR_COL,
                                         sizeof (float));
    if (hgt2 == NULL)
    {
        RETURN_ERROR ("Allocating HGT_2 memory", FUNC_NAME, FAILURE);
    }

    spfh2 = (float **) allocate_2d_array (P_LAYER, NARR_ROW * NARR_COL,
                                          sizeof (float));
    if (spfh2 == NULL)
    {
        RETURN_ERROR ("Allocating SPFH_2 memory", FUNC_NAME, FAILURE);
    }

    tmp2 = (float **) allocate_2d_array (P_LAYER, NARR_ROW * NARR_COL,
                                         sizeof (float));
    if (tmp2 == NULL)
    {
        RETURN_ERROR ("Allocating TMP_2 memory", FUNC_NAME, FAILURE);
    }

    /* Read in NARR height for time before landsat acqusition */
    if (read_narr_parameter_values (layers, "HGT_1", hgt1) != SUCCESS)
    {
        RETURN_ERROR ("Failed loading HGT_1 parameters", FUNC_NAME, FAILURE);
    }

    /* Read in NARR specific humidity for time before landsat acqusition */
    if (read_narr_parameter_values (layers, "SPFH_1", spfh1) != SUCCESS)
    {
        RETURN_ERROR ("Failed loading SPFH_1 parameters", FUNC_NAME, FAILURE);
    }

    /* Read in NARR temperature for time before landsat acqusition */
    if (read_narr_parameter_values (layers, "TMP_1", tmp1) != SUCCESS)
    {
        RETURN_ERROR ("Failed loading TMP_1 parameters", FUNC_NAME, FAILURE);
    }

    /* Read in NARR height for time after landsat acqusition */
    if (read_narr_parameter_values (layers, "HGT_2", hgt2) != SUCCESS)
    {
        RETURN_ERROR ("Failed loading HGT_2 parameters", FUNC_NAME, FAILURE);
    }

    /* Read in NARR specific humidity for time after landsat acqusition */
    if (read_narr_parameter_values (layers, "SPFH_2", spfh2) != SUCCESS)
    {
        RETURN_ERROR ("Failed loading SPFH_2 parameters", FUNC_NAME, FAILURE);
    }

    /* Read in NARR temperature for time after landsat acqusition */
    if (read_narr_parameter_values (layers, "TMP_2", tmp2) != SUCCESS)
    {
        RETURN_ERROR ("Failed loading TMP_2 parameters", FUNC_NAME, FAILURE);
    }

    /* ==================================================================== */

    /* Allocate memory for height of NARR points within the rectangular */
    narr_lat = (float *) malloc (num_points * sizeof (float));
    if (narr_lat == NULL)
    {
        RETURN_ERROR ("Allocating narr_lat memory", FUNC_NAME, FAILURE);
    }

    narr_lon = (float *) malloc (num_points * sizeof (float));
    if (narr_lon == NULL)
    {
        RETURN_ERROR ("Allocating narr_lon memory", FUNC_NAME, FAILURE);
    }

    narr_hgt1 = (float **) allocate_2d_array (P_LAYER, num_points,
                                              sizeof (float));
    if (narr_hgt1 == NULL)
    {
        RETURN_ERROR ("Allocating narr_hgt1 memory", FUNC_NAME, FAILURE);
    }

    narr_hgt2 = (float **) allocate_2d_array (P_LAYER, num_points,
                                              sizeof (float));
    if (narr_hgt2 == NULL)
    {
        RETURN_ERROR ("Allocating narr_hgt2 memory", FUNC_NAME, FAILURE);
    }

    /* Allocate memory for humidity of NARR points within the rectangular */
    narr_spfh1 = (float **) allocate_2d_array (P_LAYER, num_points,
                                               sizeof (float));
    if (narr_spfh1 == NULL)
    {
        RETURN_ERROR ("Allocating narr_spfh1 memory", FUNC_NAME, FAILURE);
    }

    narr_spfh2 = (float **) allocate_2d_array (P_LAYER, num_points,
                                               sizeof (float));
    if (narr_spfh2 == NULL)
    {
        RETURN_ERROR ("Allocating narr_spfh2 memory", FUNC_NAME, FAILURE);
    }

    /* Allocate memory for temperature of NARR points within the rectangular */
    narr_tmp1 = (float **) allocate_2d_array (P_LAYER, num_points,
                                              sizeof (float));
    if (narr_tmp1 == NULL)
    {
        RETURN_ERROR ("Allocating narr_tmp1 memory", FUNC_NAME, FAILURE);
    }

    narr_tmp2 = (float **) allocate_2d_array (P_LAYER, num_points,
                                              sizeof (float));
    if (narr_tmp2 == NULL)
    {
        RETURN_ERROR ("Allocating narr_tmp2 memory", FUNC_NAME, FAILURE);
    }

    /* extract coordinates within the NARR rectangle */
    for (j = min_jay; j <= max_jay; j++)
    {
        for (i = min_eye; i <= max_eye; i++)
        {
            narr_lat[(j - min_jay) * num_eyes + (i - min_eye)] = lat[i][j];
            narr_lon[(j - min_jay) * num_eyes + (i - min_eye)] = lon[i][j];
        }
    }
    for (k = 0; k < P_LAYER; k++)
    {
        for (j = min_jay; j <= max_jay; j++)
        {
            for (i = min_eye; i <= max_eye; i++)
            {
                narr_hgt1[k][(j - min_jay) * num_eyes + (i - min_eye)]
                    = hgt1[k][j * NARR_ROW + i];
                narr_spfh1[k][(j - min_jay) * num_eyes + (i - min_eye)]
                    = spfh1[k][j * NARR_ROW + i];
                narr_tmp1[k][(j - min_jay) * num_eyes + (i - min_eye)]
                    = tmp1[k][j * NARR_ROW + i];
                narr_hgt2[k][(j - min_jay) * num_eyes + (i - min_eye)]
                    = hgt2[k][j * NARR_ROW + i];
                narr_spfh2[k][(j - min_jay) * num_eyes + (i - min_eye)]
                    = spfh2[k][j * NARR_ROW + i];
                narr_tmp2[k][(j - min_jay) * num_eyes + (i - min_eye)]
                    = tmp2[k][j * NARR_ROW + i];
            }
        }
    }

    /* Release memory */
    if (free_2d_array ((void **) lat) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: lat\n", FUNC_NAME, FAILURE);
    }
    lat = NULL;

    if (free_2d_array ((void **) lon) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: lon\n", FUNC_NAME, FAILURE);
    }
    lon = NULL;

    if (free_2d_array ((void **) hgt1) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: hgt1\n", FUNC_NAME, FAILURE);
    }
    hgt1 = NULL;

    if (free_2d_array ((void **) spfh1) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: spfh1\n", FUNC_NAME, FAILURE);
    }
    spfh1 = NULL;

    if (free_2d_array ((void **) tmp1) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: tmp1\n", FUNC_NAME, FAILURE);
    }
    tmp1 = NULL;

    if (free_2d_array ((void **) hgt2) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: ght2\n", FUNC_NAME, FAILURE);
    }
    hgt2 = NULL;

    if (free_2d_array ((void **) spfh2) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: spfh2\n", FUNC_NAME, FAILURE);
    }
    spfh2 = NULL;

    if (free_2d_array ((void **) tmp2) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: tmp2\n", FUNC_NAME, FAILURE);
    }
    tmp2 = NULL;

    /* ==================================================================== */

    /* Allocate memory */
    pressure = (float **) allocate_2d_array (P_LAYER, num_points,
                                             sizeof (float));
    if (pressure == NULL)
    {
        RETURN_ERROR ("Allocating pressure memory", FUNC_NAME, FAILURE);
    }

    for (i = 0; i < P_LAYER; i++)
    {
        for (j = 0; j < num_points; j++)
        {
            pressure[i][j] = layers[i];
        }
    }

    /* ==================================================================== */

    narr_height1 = (float **) allocate_2d_array (P_LAYER, num_points,
                                                 sizeof (float));
    if (narr_height1 == NULL)
    {
        RETURN_ERROR ("Allocating narr_height1 memory", FUNC_NAME, FAILURE);
    }

    narr_height2 = (float **) allocate_2d_array (P_LAYER, num_points,
                                                 sizeof (float));
    if (narr_height2 == NULL)
    {
        RETURN_ERROR ("Allocating narr_height2 memory", FUNC_NAME, FAILURE);
    }

    narr_rh1 = (float **) allocate_2d_array (P_LAYER, num_points,
                                             sizeof (float));
    if (narr_rh1 == NULL)
    {
        RETURN_ERROR ("Allocating narr_rh1 memory", FUNC_NAME, FAILURE);
    }

    narr_rh2 = (float **) allocate_2d_array (P_LAYER, num_points,
                                             sizeof (float));
    if (narr_rh2 == NULL)
    {
        RETURN_ERROR ("Allocating narr_rh2 memory", FUNC_NAME, FAILURE);
    }

    /* convert grib data to variables to be input to MODTRAN */
    if (convert_geopotential_geometric (num_points, narr_lat,
                                        narr_hgt1, narr_height1) != SUCCESS)
    {
        RETURN_ERROR ("Calling convert_geopotential_geometric for height 1",
                      FUNC_NAME, FAILURE);
    }

    if (convert_geopotential_geometric (num_points, narr_lat,
                                        narr_hgt2, narr_height2) != SUCCESS)
    {
        RETURN_ERROR ("Calling convert_geopotential_geometric for height 2",
                      FUNC_NAME, FAILURE);
    }

    if (convert_sh_rh (num_points, narr_lat, narr_spfh1, narr_tmp1,
                       pressure, narr_rh1) != SUCCESS)
    {
        RETURN_ERROR ("Calling convert_sh_rh for temp 1", FUNC_NAME,
                      FAILURE);
    }

    if (convert_sh_rh (num_points, narr_lat, narr_spfh2, narr_tmp2,
                       pressure, narr_rh2) != SUCCESS)
    {
        RETURN_ERROR ("Calling convert_sh_rh for temp 2", FUNC_NAME,
                      FAILURE);
    }

    /* free allocated memory */
    if (free_2d_array ((void **) narr_hgt1) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: narr_hgt1\n", FUNC_NAME, FAILURE);
    }
    narr_hgt1 = NULL;

    if (free_2d_array ((void **) narr_hgt2) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: narr_hgt2\n", FUNC_NAME, FAILURE);
    }
    narr_hgt2 = NULL;

    if (free_2d_array ((void **) narr_spfh1) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: narr_spfh1\n", FUNC_NAME, FAILURE);
    }
    narr_spfh1 = NULL;

    if (free_2d_array ((void **) narr_spfh2) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: narr_spfh2\n", FUNC_NAME, FAILURE);
    }
    narr_spfh2 = NULL;

    /* ==================================================================== */

    /* determine three hour-increment before and after scene center scan time */
    /* TODO TODO TODO - Does not take into consideration day traversal..... */
    /* TODO TODO TODO - Does not take into consideration day traversal..... */
    /* TODO TODO TODO - Does not take into consideration day traversal..... */
    rem1 = input->meta.acq_date.hour % 3;
    rem2 = 3 - rem1;
    hour1 = (float) (input->meta.acq_date.hour - rem1);
    hour2 = (float) (input->meta.acq_date.hour + rem2);
    inv_hour_diff = 1.0 / (hour2 - hour1);

    /* Round to the nearest minute */
    if ((input->meta.acq_date.second - 30.0) >= MINSIGMA)
        input->meta.acq_date.minute++;

    /* convert hour-min acquisition time to decimal time */
    time = (float) input->meta.acq_date.hour
           + (float) input->meta.acq_date.minute / 60.0F;
    time_diff = time - hour1;

    /* ==================================================================== */

    /* Allocate memory */
    narr_height = (float **) allocate_2d_array (P_LAYER, num_points,
                                                sizeof (float));
    if (narr_height == NULL)
    {
        RETURN_ERROR ("Allocating narr_height memory", FUNC_NAME, FAILURE);
    }

    /* Allocate memory */
    narr_rh = (float **) allocate_2d_array (P_LAYER, num_points,
                                            sizeof (float));
    if (narr_rh == NULL)
    {
        RETURN_ERROR ("Allocating narr_rh memory", FUNC_NAME, FAILURE);
    }

    /* Allocate memory */
    narr_tmp = (float **) allocate_2d_array (P_LAYER, num_points,
                                             sizeof (float));
    if (narr_tmp == NULL)
    {
        RETURN_ERROR ("Allocating narr_tmp memory", FUNC_NAME, FAILURE);
    }

    /* linearly interpolate geometric height, relative humidity, and
       temperature for NARR points within Landsat scene this is the NARR data
       corresponding to the acquisition time of the Landsat image converted to
       appropriated variable for MODTRAN input */
    for (i = 0; i < P_LAYER; i++)
    {
        for (j = 0; j < num_points; j++)
        {
            narr_height[i][j] = narr_height1[i][j] + time_diff
                                * ((narr_height2[i][j] - narr_height1[i][j])
                                   * inv_hour_diff);
            narr_rh[i][j] = narr_rh1[i][j] + time_diff
                            * ((narr_rh2[i][j] - narr_rh1[i][j])
                               * inv_hour_diff);
            narr_tmp[i][j] = narr_tmp1[i][j] + time_diff
                             * ((narr_tmp2[i][j] - narr_tmp1[i][j])
                                * inv_hour_diff);
        }
    }

    /* Free allocated memory */
    if (free_2d_array ((void **) narr_height1) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: narr_height1\n", FUNC_NAME, FAILURE);
    }
    narr_height1 = NULL;

    if (free_2d_array ((void **) narr_height2) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: narr_height2\n", FUNC_NAME, FAILURE);
    }
    narr_height2 = NULL;

    if (free_2d_array ((void **) narr_rh1) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: narr_rh1\n", FUNC_NAME, FAILURE);
    }
    narr_rh1 = NULL;

    if (free_2d_array ((void **) narr_rh2) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: narr_rh2\n", FUNC_NAME, FAILURE);
    }
    narr_rh2 = NULL;

    if (free_2d_array ((void **) narr_tmp1) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: narr_tmp1\n", FUNC_NAME, FAILURE);
    }
    narr_tmp1 = NULL;

    if (free_2d_array ((void **) narr_tmp2) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: narr_tmp2\n", FUNC_NAME, FAILURE);
    }
    narr_tmp2 = NULL;

    /************************** Build tape 5 files **************************/

    /* Allocate memory */
    stan_height = (float *) malloc (STAN_LAYER * sizeof (float));
    if (stan_height == NULL)
    {
        RETURN_ERROR ("Allocating stan_height memory", FUNC_NAME, FAILURE);
    }

    stan_pre = (float *) malloc (STAN_LAYER * sizeof (float));
    if (stan_pre == NULL)
    {
        RETURN_ERROR ("Allocating stan_pre memory", FUNC_NAME, FAILURE);
    }

    stan_temp = (float *) malloc (STAN_LAYER * sizeof (float));
    if (stan_temp == NULL)
    {
        RETURN_ERROR ("Allocating stan_temp memory", FUNC_NAME, FAILURE);
    }

    stan_rh = (float *) malloc (STAN_LAYER * sizeof (float));
    if (stan_rh == NULL)
    {
        RETURN_ERROR ("Allocating stan_rh memory", FUNC_NAME, FAILURE);
    }

    if (read_std_mid_lat_summer_atmos (lst_data_dir, stan_height, stan_pre,
                                       stan_temp, stan_rh)
        != SUCCESS)
    {
        RETURN_ERROR ("Failed loading std_mid_lat_summer_atmos",
                      FUNC_NAME, FAILURE);
    }

    /* determine number of MODTRAN runs */
    num_modtran_runs = num_points * NUM_ELEVATIONS * 3;

    /* Allocate memory */
    counter = (int *) malloc (STAN_LAYER * sizeof (int));
    if (counter == NULL)
    {
        RETURN_ERROR ("Allocating counter memory", FUNC_NAME, FAILURE);
    }

    case_list = (char **) allocate_2d_array (num_modtran_runs, MAX_STR_LEN,
                                             sizeof (char));
    if (case_list == NULL)
    {
        RETURN_ERROR ("Allocating case_list memory", FUNC_NAME, FAILURE);
    }

    command_list = (char **) allocate_2d_array (num_modtran_runs, MAX_STR_LEN,
                                                sizeof (char));
    if (command_list == NULL)
    {
        RETURN_ERROR ("Allocating command_list memory", FUNC_NAME, FAILURE);
    }

    modtran_path = getenv ("MODTRAN_PATH");
    if (modtran_path == NULL)
    {
        RETURN_ERROR ("MODTRAN_PATH environment variable is not set",
                      FUNC_NAME, FAILURE);
    }

    modtran_data_dir = getenv ("MODTRAN_DATA_DIR");
    if (modtran_data_dir == NULL)
    {
        RETURN_ERROR ("MODTRAN_DATA_DIR environment variable is not set",
                      FUNC_NAME, FAILURE);
    }

    getcwd (curr_path, MAX_STR_LEN);

    /* Allocate some temp memory */
    temp_height = (float *) malloc (MAX_MODTRAN_LAYER * sizeof (float));
    if (temp_height == NULL)
    {
        RETURN_ERROR ("Allocating temp_height memory", FUNC_NAME, FAILURE);
    }

    temp_pressure = (float *) malloc (MAX_MODTRAN_LAYER * sizeof (float));
    if (temp_pressure == NULL)
    {
        RETURN_ERROR ("Allocating temp_pressure memory", FUNC_NAME, FAILURE);
    }

    temp_temp = (float *) malloc (MAX_MODTRAN_LAYER * sizeof (float));
    if (temp_temp == NULL)
    {
        RETURN_ERROR ("Allocating temp_temp memory", FUNC_NAME, FAILURE);
    }

    temp_rh = (float *) malloc (MAX_MODTRAN_LAYER * sizeof (float));
    if (temp_rh == NULL)
    {
        RETURN_ERROR ("Allocating temp_rh memory", FUNC_NAME, FAILURE);
    }

    for (i = 0; i < num_points; i++)
    {
        /* ****************************************************************
           MODTRAN uses longitudinal degree values from 0 to 360 starting at
           Greenwich and moving west.

           The following logic fixes the longitude to be for MODTRAN and we
           also use those values for generating the output directory names.
           **************************************************************** */
        if (narr_lon[i] < 0)
        {
            /* We are a west longitude value so negate it to the positive
               equivalent value.
               i.e.   W40 normally represented as -40 would be changed to a
                      positive 40 value. */
            narr_lon[i] = -narr_lon[i];
        }
        else
        {
            /* We are a east longitude value so fix it to be greater than
               180 west.
               i.e.   E10 would be turned into W350, and be positive 350 not
                      negative.  */
            narr_lon[i] = 360.0 - narr_lon[i];
        }

        /* Figure out the lat and lon strings to use */
        if ((fabs (narr_lat[i]) - 100.0) < MINSIGMA)
            snprintf (lat_str, sizeof (lat_str), "%6.3f", narr_lat[i]);
        else
            snprintf (lat_str, sizeof (lat_str), "%6.2f", narr_lat[i]);

        if ((fabs (narr_lon[i]) - 100.0) < MINSIGMA)
            snprintf (lon_str, sizeof (lon_str), "%6.3f", narr_lon[i]);
        else
            snprintf (lon_str, sizeof (lon_str), "%6.2f", narr_lon[i]);

        /* Create the name of the directory for the current NARR point */
        snprintf (current_point, sizeof (current_point),
                  "%s_%s", lat_str, lon_str);

        /* Create the directory */
        snprintf(msg_str, sizeof (msg_str),
                 "Creating directory [%s]", current_point);
        LOG_MESSAGE(msg_str, FUNC_NAME);
        if (mkdir (current_point, 0755) != SUCCESS)
        {
            if (errno != EEXIST)
            {
                snprintf (msg_str, sizeof (msg_str),
                          "Failed creating directory [%s]", current_point);
                RETURN_ERROR (msg_str, FUNC_NAME, FAILURE);
            }
        }

        /* determine latitude and longitude of current NARR point and insert
           into tail file */
        /* TODO - The sed commands could be combined into one to save
                  traversing the file three times.  But the file is small. */
        snprintf (command, sizeof (command),
                  "cat %s/modtran_tail.txt"
                  " | sed 's/latitu/%s/'"
                  " | sed 's/longit/%s/'"
                  " | sed 's/jay/%d/' > newTail.txt",
                  lst_data_dir, lat_str, lon_str, input->meta.acq_date.doy);
        if (system (command) != SUCCESS)
        {
            RETURN_ERROR ("Failed creating newTail.txt", FUNC_NAME, FAILURE);
        }

        /* Clear the temp memory */
        for (j = 0; j < MAX_MODTRAN_LAYER; j++)
        {
            temp_height[j] = 0.0;
            temp_pressure[j] = 0.0;
            temp_temp[j] = 0.0;
            temp_rh[j] = 0.0;
        }

        /* set lowest altitude is the first geometric height at that NARR
           point (if positive) and (if negative set to zero) */
        if (narr_height[0][i] < 0)
            gndalt[0] = 0.0;
        else
            gndalt[0] = narr_height[0][i];

        /* iterate through all ground altitudes at which MODTRAN is run */
        for (j = 0; j < NUM_ELEVATIONS; j++)
        {
            /* create a directory for the current height */
            snprintf (current_gdalt, sizeof (current_gdalt),
                      "%s/%5.3f", current_point, gndalt[j]);

            /* Create the directory */
            if (mkdir (current_gdalt, 0755) != SUCCESS)
            {
                if (errno != EEXIST)
                {
                    snprintf (msg_str, sizeof (msg_str),
                              "Failed creating directory [%s]", current_gdalt);
                    RETURN_ERROR (msg_str, FUNC_NAME, FAILURE);
                }
            }

            /* determine layers below current gndalt and closest index
               above and below */
            for (k = 0; k < P_LAYER; k++)
            {
                if ((narr_height[k][i] - gndalt[j]) >= MINSIGMA)
                {
                    index_below = k - 1;
                    index_above = k;
                    break;
                }
            }

            if (index_below < 0)
            {
                index_below = 0;
                index_above = 1;
            }

            /* To save divisions */
            inv_narr_height_diff = 1.0 / (narr_height[index_above][i]
                                          - narr_height[index_below][i]);

            /* linearly interpolate pressure, temperature, and relative 
               humidity to gndalt for lowest layer */
            new_pressure = pressure[index_below][i]
                           + (gndalt[j] - narr_height [index_below][i])
                           * ((pressure[index_above][i]
                               - pressure[index_below][i])
                              * inv_narr_height_diff);
            new_temp = narr_tmp[index_below][i]
                       + (gndalt[j] - narr_height[index_below][i])
                       * ((narr_tmp[index_above][i] - narr_tmp[index_below][i])
                          * inv_narr_height_diff);
            new_rh = narr_rh[index_below][i]
                     + (gndalt[j] - narr_height[index_below][i])
                     * ((narr_rh[index_above][i] - narr_rh[index_below][i])
                        * inv_narr_height_diff);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /* create arrays containing only layers to be included in current
               tape5 file */
            index = 0;
            temp_height[index] = gndalt[j];
            temp_pressure[index] = new_pressure;
            temp_temp[index] = new_temp;
            temp_rh[index] = new_rh;
            index++;

            for (k = index_above; k < P_LAYER; k++)
            {
                temp_height[index] = narr_height[k][i];
                temp_pressure[index] = pressure[k][i];
                temp_temp[index] = narr_tmp[k][i];
                temp_rh[index] = narr_rh[k][i];
                index++;
            }


            /* MODTRAN throws an error when there are two identical layers in
               the tape5 file, if the current ground altitude and the next
               highest layer are close enough, eliminate interpolated layer */
            if (fabs (gndalt[j] - narr_height[index_above][i] - 0.001)
                < MINSIGMA)
            {
                index = 0;
                for (k = index_above; k < P_LAYER; k++)
                {
                    temp_height[index] = narr_height[k][i];
                    temp_pressure[index] = pressure[k][i];
                    temp_temp[index] = narr_tmp[k][i];
                    temp_rh[index] = narr_rh[k][i];
                    index++;
                }
            }

            /* determine maximum height of NARR layers and where the standard 
               atmosphere is greater than this */
            index2 = 0;
            for (k = 0; k < STAN_LAYER; k++)
            {
                if (stan_height[k] > narr_height[P_LAYER - 1][i])
                {
                    counter[index2] = k;
                    index2++;
                }
            }

            /* if there are at least three layers above to highest NARR layer,
               add standard atmosphere layers and linearly interpolate height,
               pressure, temp, and rel hum to create a smooth transition
               between the NARR layers and the standard upper atmosphere */
// RDD - Did this to remove a compiler warning, but something else may need to be done.
            new_height = 0;
            if (index2 >= 3)
            {
                new_height = (stan_height[counter[2]]
                              + temp_height[index - 1]) / 2.0;
                new_pressure = temp_pressure[index - 1]
                               + (new_height - temp_height[index - 1])
                               * ((stan_pre[counter[2]]
                                   - temp_pressure[index - 1])
                                  / (stan_height[counter[2]]
                                     - temp_height[index - 1]));
                new_temp = temp_temp[index - 1]
                           + (new_height - temp_height[index - 1])
                           * ((stan_temp[counter[2]] - temp_temp[index - 1])
                              / (stan_height[counter[2]]
                                 - temp_height[index - 1]));
                new_rh = temp_rh[counter[2]]
                         + (new_height - temp_height[index - 1])
                         * ((stan_rh[counter[2]] - temp_rh[index - 1])
                            / (stan_height[counter[2]]
                               - temp_height[index - 1]));
            }

            /* concatenate NARR layers, new layer, and standard atmosphere
               layers */
            temp_height[index] = new_height;
            temp_pressure[index] = new_pressure;
            temp_temp[index] = new_temp;
            temp_rh[index] = new_rh;
            index++;
            for (k = 2; k < index2; k++)
            {
                temp_height[index] = stan_height[counter[k]];
                temp_pressure[index] = stan_pre[counter[k]];
                temp_temp[index] = stan_temp[counter[k]];
                temp_rh[index] = stan_rh[counter[k]];
                index++;
            }

            /* write atmospheric layers to a text file in format proper for
               tape5 file */
            fd = fopen ("tempLayers.txt", "w");
            if (fd == NULL)
            {
                RETURN_ERROR ("Opening file: tempLayers.txt\n",
                              FUNC_NAME, FAILURE);
            }

            /* Write out the intermediate file */
            for (k = 0; k < index; k++)
            {
                fprintf (fd, "%10.3f%10.3e%10.3e%10.3e%10.3e%10.3e%16s\n",
                         temp_height[k], temp_pressure[k], temp_temp[k],
                         temp_rh[k], 0.0, 0.0, "AAH             ");
            }

            /* Close the intermediate file */
            if (fclose (fd) != SUCCESS)
            {
                RETURN_ERROR ("Closing file: tempLayers.txt\n",
                              FUNC_NAME, FAILURE);
            }

            /* determine number of layers for current ground altitude and
               insert into head file */
            sprintf (command, "cat %s/modtran_head.txt"
                              " | sed 's/nml/%d/'" " > newHead.txt",
                              lst_data_dir, index);

            if (system (command) != SUCCESS)
            {
                RETURN_ERROR ("system call 5", FUNC_NAME, FAILURE);
            }

            /* insert current ground altitude into head file */
            sprintf (command, "cat newHead.txt"
                     " | sed 's/gdalt/%5.3f/'" " > newHead2.txt", gndalt[j]);
            if (system (command) != SUCCESS)
            {
                RETURN_ERROR ("system call 6", FUNC_NAME, FAILURE);
            }

            /* iterate through [temperature,albedo] pairs at which to run
               MODTRAN */
            for (k = 0; k <= 2; k++)
            {
                if (k == 2)
                    sprintf (temp_out, "%s", "000");
                else
                    sprintf (temp_out, "%d", (int) tmp[k]);

                /* create directory for the current temperature */
                sprintf (current_temp, "%s/%s", current_gdalt, temp_out);

                if (mkdir (current_temp, 0755) != SUCCESS)
                {
                    if (errno != EEXIST)
                    {
                        snprintf (msg_str, sizeof (msg_str),
                                  "Failed creating directory [%s]",
                                  current_temp);
                        RETURN_ERROR (msg_str, FUNC_NAME, FAILURE);
                    }
                }

                /* insert current temperature into head file */
                sprintf (command, "cat newHead2.txt"
                         " | sed 's/tmp/%s/'" " > newHead3.txt", temp_out);
                if (system (command) != SUCCESS)
                {
                    RETURN_ERROR ("system call 8", FUNC_NAME, FAILURE);
                }

                /* create directory for the current albedo */
                sprintf (current_alb, "%s/%3.1f", current_temp, alb[k]);
                if (mkdir (current_alb, 0755) != SUCCESS)
                {
                    if (errno != EEXIST)
                    {
                        snprintf (msg_str, sizeof (msg_str),
                                  "Failed creating directory [%s]",
                                  current_alb);
                        RETURN_ERROR (msg_str, FUNC_NAME, FAILURE);
                    }
                }

                /* insert current albedo into head file */
                sprintf (command, "cat newHead3.txt"
                         " | sed 's/alb/%4.2f/'" " > newHead4.txt", alb[k]);
                if (system (command) != SUCCESS)
                {
                    RETURN_ERROR ("system call 10", FUNC_NAME, FAILURE);
                }

                /* concatenate head file, atmospheric layers, and tail file to
                   create a tape5 file for MODTRAN specific to this location
                   and ground altitude with variables for temperature and
                   albedo */
                sprintf (command, "cat newHead4.txt"
                         " tempLayers.txt"
                         " newTail.txt" " > %s/tape5", current_alb);
                if (system (command) != SUCCESS)
                {
                    RETURN_ERROR ("system call 11", FUNC_NAME, FAILURE);
                }

                /* create string for case list containing location of current
                   tape5 file
                   create string for command list containing commands for
                   MODTRAN run
                   iterate entry count */
                case_counter = i * NUM_ELEVATIONS * 3 + j * 3 + k;
                sprintf (case_list[case_counter], "%s/%s",
                         curr_path, current_alb);
                sprintf (command_list[case_counter],
                         "pushd %s; ln -s %s; %s/Mod90_5.2.2.exe; popd",
                         case_list[case_counter], modtran_data_dir,
                         modtran_path);
            } /* END - Tempurature Albedo Pairs */
        } /* END - ground altitude ran by MODTRAN */
    } /* END - number of points */

    /* Free the temp memory */
    free(temp_height);
    free(temp_pressure);
    free(temp_temp);
    free(temp_rh);

    /* write caseList to a file */
    fd = fopen ("caseList", "w");
    if (fd == NULL)
    {
        RETURN_ERROR ("Opening file: caseList\n", FUNC_NAME, FAILURE);
    }

    /* Write out the caseList file */
    for (k = 0; k < num_points * NUM_ELEVATIONS * 3; k++)
    {
        fprintf (fd, "%s\n", case_list[k]);
    }

    /* Close the caseList file */
    if (fclose (fd) != SUCCESS)
    {
        RETURN_ERROR ("Closing file: caseList\n", FUNC_NAME, FAILURE);
    }

    /* write commandList to a file */
    fd = fopen ("commandList", "w");
    if (fd == NULL)
    {
        RETURN_ERROR ("Opening file: commandList\n", FUNC_NAME, FAILURE);
    }

    /* Write out the commandList file */
    for (k = 0; k < num_points * NUM_ELEVATIONS * 3; k++)
    {
        fprintf (fd, "%s\n", command_list[k]);
    }

    /* Close the commandList file */
    if (fclose (fd) != SUCCESS)
    {
        RETURN_ERROR ("Closing file: commandList\n", FUNC_NAME, FAILURE);
    }

    /* Free memory allocation */
    if (free_2d_array ((void **) pressure) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: pressure\n", FUNC_NAME, FAILURE);
    }

    if (free_2d_array ((void **) narr_height) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: narr_height\n", FUNC_NAME, FAILURE);
    }

    if (free_2d_array ((void **) narr_rh) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: narr_rh\n", FUNC_NAME, FAILURE);
    }

    if (free_2d_array ((void **) narr_tmp) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: narr_tmp\n", FUNC_NAME, FAILURE);
    }

    /* Free memory allocation */
    if (free_2d_array ((void **) case_list) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: case_list\n", FUNC_NAME, FAILURE);
    }

    if (free_2d_array ((void **) command_list) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: command_list\n", FUNC_NAME, FAILURE);
    }

    /* Assign to the output variable */
    *num_pts = num_points;

    return SUCCESS;
}
