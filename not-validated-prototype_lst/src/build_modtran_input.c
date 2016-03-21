
#include <unistd.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <sys/stat.h>
#include <errno.h>


#include "const.h"
#include "2d_array.h"
#include "utilities.h"
#include "input.h"
#include "lst_types.h"
#include "build_points.h"


#define STANRDARD_GRAVITY_IN_M_PER_SEC_SQRD 9.80665
#define EQUATORIAL_RADIUS_IN_KM (UTM_EQUATORIAL_RADIUS / 1000.0)
#define POLAR_RADIUS_IN_KM (UTM_POLAR_RADIUS / 1000.0)


#define P_LAYER 29
#define STANDARD_LAYERS 30
#define MAX_MODTRAN_LAYER 150


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
    int num_points,         /* I: number of points */
    double *lat,            /* I: latitude in degrees */
    double **geo_potential, /* I: geo_potential height */
    double **geo_metric     /* O: geo_metric height */
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

    /* Calculate the geometric height for each point in each layer */
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
    double *lat,       /* I: latitude in degrees */
    double **spec_hum,
    double **temp_k,
    double **pressure,
    double **rh        /* O: relative humidity */
)
{
    char FUNC_NAME[] = "convert_sh_rh";
    int layer;
    int point;
    double mh20 = 18.01534;
    double mdry = 28.9644;

    double **goff;
    double **ph20;

    /* Allocate memory */
    goff = (double **) allocate_2d_array (P_LAYER, num_points, sizeof (double));
    if (goff == NULL)
    {
        RETURN_ERROR ("Allocating goff memory", FUNC_NAME, FAILURE);
    }

    ph20 = (double **) allocate_2d_array (P_LAYER, num_points, sizeof (double));
    if (ph20 == NULL)
    {
        RETURN_ERROR ("Allocating  memory", FUNC_NAME, FAILURE);
    }

    for (layer = 0; layer < P_LAYER; layer++)
    {
        for (point = 0; point < num_points; point++)
        {
            /* calculate vapor pressure at given temperature - hpa */
            goff[layer][point] = -7.90298 * (373.16 / temp_k[layer][point]
                                             - 1.0)
                                 + 5.02808 * log10 (373.16
                                                    / temp_k[layer][point])
                                 - 1.3816e-7
                                 * (pow (10.0, (11.344
                                                * (1.0 - (temp_k[layer][point]
                                                          / 373.16))))
                                    - 1.0)
                                 + 8.1328e-3
                                 * (pow (10.0, (-3.49149
                                               * (373.16 / temp_k[layer][point]
                                                  - 1.0)))
                                        - 1.0)
                                 + log10 (1013.246); /* hPa */

            /* calculate partial pressure */
            ph20[layer][point] = (spec_hum[layer][point]
                                  * pressure[layer][point] * mdry)
                                 / (mh20
                                    - spec_hum[layer][point] * mh20
                                    + spec_hum[layer][point] * mdry);

            /* calculate relative humidity */
            rh[layer][point] = (ph20[layer][point]
                                / pow (10.0, goff[layer][point])) * 100.0;
        }
    }

    /* Free allocated memory */
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
MODULE:  read_std_mid_lat_summer_atmos

PURPOSE: Reads the standard atmosphere into memory

RETURN: SUCCESS
        FAILURE
******************************************************************************/
int read_std_mid_lat_summer_atmos
(
    char *lst_data_dir,
    double *stan_height,
    double *stan_pre,
    double *stan_temp,
    double *stan_rh
)
{
    char FUNC_NAME[] = "read_std_mid_lat_summer_atmos";
    int layer;
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

    for (layer = 0; layer < STANDARD_LAYERS; layer++)
    {
        if (fscanf (fd, "%lf %lf %lf %lf",
                    &stan_height[layer], &stan_pre[layer],
                    &stan_temp[layer], &stan_rh[layer]) == EOF)
        {
            RETURN_ERROR ("End of file (EOF) is met before STANDARD_LAYERS"
                          " lines", FUNC_NAME, FAILURE);
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

PURPOSE: Reads the NARR parameter values into memory

RETURN: SUCCESS
        FAILURE
******************************************************************************/
int read_narr_parameter_values
(
    int *layers,
    char *parameter,
    double **output_2d_array
)
{
    char FUNC_NAME[] = "read_narr_parameter_values";
    char msg_str[MAX_STR_LEN];
    char parm_filename[PATH_MAX];
    int layer;
    int point;
    int count;
    int file_rows;
    int file_cols;
    FILE *fd = NULL;

    /* Read each layers parameter file into memory */
    for (layer = 0; layer < P_LAYER; layer++)
    {
        /* Build the full path to the parameter file */
        count = snprintf (parm_filename, sizeof (parm_filename),
                           "%s/%d%s", parameter, layers[layer], ".txt");
        if (count < 0 || count >= sizeof (parm_filename))
        {
            RETURN_ERROR ("Failed initializing parm_filename variable for"
                          " NARR parameters", FUNC_NAME, FAILURE);
        }

#if 0
        snprintf(msg_str, sizeof (msg_str), "Reading [%s]", parm_filename);
        LOG_MESSAGE (msg_str, FUNC_NAME);
#endif

        /* Open the parameter file for reading */
        fd = fopen (parm_filename, "r");
        if (fd == NULL)
        {
            RETURN_ERROR ("Can't HGT_1 txt file", FUNC_NAME, FAILURE);
        }

        /* Read the number of rows and columns and verify them */
        count = fscanf (fd, "%d %d", &file_cols, &file_rows);
        if (count != 2)
        {
            RETURN_ERROR ("Failed reading columns and rows",
                          FUNC_NAME, FAILURE);
        }
        if (file_rows != NARR_ROWS || file_cols != NARR_COLS)
        {
            RETURN_ERROR ("Parameter file contains invalid number of rows and"
                          " columns", FUNC_NAME, FAILURE);
        }

        /* Read the values into memory */
        for (point = 0; point < NARR_ROWS * NARR_COLS; point++)
        {
            if (fscanf (fd, "%lf", &output_2d_array[layer][point]) == EOF)
            {
                RETURN_ERROR ("End of file (EOF) is met before "
                              "NARR_ROWS * NARR_COLS lines",
                               FUNC_NAME, FAILURE);
            }
        }

        /* Close the current parameter file */
        if (fclose (fd) != SUCCESS)
        {
            snprintf(msg_str, sizeof (msg_str),
                     "Closing file: %s/%d.txt\n", parameter, layers[layer]);
            RETURN_ERROR (msg_str, FUNC_NAME, FAILURE);
        }
    }

    return SUCCESS;
}


/******************************************************************************
MODULE:  build_modtran_input

PURPOSE: Creates directories and writes tape5 files

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/21/2014   Song Guo         Original Development
******************************************************************************/
int build_modtran_input
(
    Input_Data_t *input,       /* I: input structure */
    REANALYSIS_POINTS *points, /* I/O: The coordinate points */
    bool verbose,         /* I: value to indicate if intermediate messages
                                should be printed */
    bool debug            /* I: value to indicate if debug should be
                                generated */
)
{
    char FUNC_NAME[] = "build_modtran_input";

    double **hgt1;
    double **spfh1;
    double **tmp1;
    double **hgt2;
    double **spfh2;
    double **tmp2;
    double **narr_hgt1;
    double **narr_spfh1;
    double **narr_tmp1;
    double **narr_hgt2;
    double **narr_spfh2;
    double **narr_tmp2;
    double **narr_height;
    double **narr_height1;
    double **narr_height2;
    double **narr_rh;
    double **narr_rh1;
    double **narr_rh2;
    double **narr_tmp;
    int row;
    int col;
    int layer;
    int point;
    int elevation;
    int temperature;
    int layers[P_LAYER] = { 1000, 975, 950, 925, 900,
                            875, 850, 825, 800, 775,
                            750, 725, 700, 650, 600,
                            550, 500, 450, 400, 350,
                            300, 275, 250, 225, 200,
                            175, 150, 125, 100 };
    double **pressure;
    int rem1;
    int rem2;
    double hour1;
    double hour2;
    double inv_hour_diff;
    double time;
    double time_diff;
    FILE *fd;
    FILE *point_list_fd;
    double *stan_height;
    double *stan_pre;
    double *stan_temp;
    double *stan_rh;
    double *temp_height;
    double *temp_pressure;
    double *temp_temp;
    double *temp_rh;
    double gndalt[NUM_ELEVATIONS] = { 0.0, 0.6, 1.1, 1.6, 2.1,
                                     2.6, 3.1, 3.6, 4.05 };
    double inv_height_diff;
    int num_modtran_runs;
    char command[PATH_MAX];
    char current_gdalt[PATH_MAX];
    char current_temp[PATH_MAX];
    char current_alb[PATH_MAX];
    char current_point[PATH_MAX];
    char curr_path[PATH_MAX];
    int layer_below = 0;
    int layer_above = NUM_ELEVATIONS;
    double new_height;
    double new_pressure;
    double new_temp;
    double new_rh;
    int curr_layer, std_layer;
    int counter[STANDARD_LAYERS];
    char temp_strs[3][4] = { "273\0", "310\0", "000\0" };
    double alb[3] = { 0.0, 0.0, 0.1 };
    char *lst_data_dir = NULL;
    char *modtran_path = NULL;
    char *modtran_data_dir = NULL;
    int case_counter;
    char lat_str[7]; /* 6 plus the string termination character */
    char lon_str[7]; /* 6 plus the string termination character */
    char msg_str[MAX_STR_LEN];


    /* Use local variables for cleaner code */
    int min_row = points->min_row;
    int max_row = points->max_row;
    int min_col = points->min_col;
    int max_col = points->max_col;
    int num_cols = points->num_cols;
    int num_points = points->num_points;

    /* Grab the environment path to the LST_DATA_DIR */
    lst_data_dir = getenv ("LST_DATA_DIR");
    if (lst_data_dir == NULL)
    {
        RETURN_ERROR ("LST_DATA_DIR environment variable is not set",
                      FUNC_NAME, FAILURE);
    }

    /* ==================================================================== */

    /* Dynamic allocate the 2d memory for time before Landsat Acq. */
    hgt1 = (double **) allocate_2d_array (P_LAYER, NARR_ROWS * NARR_COLS,
                                         sizeof (double));
    if (hgt1 == NULL)
    {
        RETURN_ERROR ("Allocating HGT_1 memory", FUNC_NAME, FAILURE);
    }

    spfh1 = (double **) allocate_2d_array (P_LAYER, NARR_ROWS * NARR_COLS,
                                          sizeof (double));
    if (spfh1 == NULL)
    {
        RETURN_ERROR ("Allocating SPFH_1 memory", FUNC_NAME, FAILURE);
    }

    tmp1 = (double **) allocate_2d_array (P_LAYER, NARR_ROWS * NARR_COLS,
                                         sizeof (double));
    if (tmp1 == NULL)
    {
        RETURN_ERROR ("Allocating TMP_1 memory", FUNC_NAME, FAILURE);
    }

    /* Dynamic allocate the 2d memory for time after Landsat Acq. */
    hgt2 = (double **) allocate_2d_array (P_LAYER, NARR_ROWS * NARR_COLS,
                                         sizeof (double));
    if (hgt2 == NULL)
    {
        RETURN_ERROR ("Allocating HGT_2 memory", FUNC_NAME, FAILURE);
    }

    spfh2 = (double **) allocate_2d_array (P_LAYER, NARR_ROWS * NARR_COLS,
                                          sizeof (double));
    if (spfh2 == NULL)
    {
        RETURN_ERROR ("Allocating SPFH_2 memory", FUNC_NAME, FAILURE);
    }

    tmp2 = (double **) allocate_2d_array (P_LAYER, NARR_ROWS * NARR_COLS,
                                         sizeof (double));
    if (tmp2 == NULL)
    {
        RETURN_ERROR ("Allocating TMP_2 memory", FUNC_NAME, FAILURE);
    }

    /* Read in NARR height for time before Landsat acqusition */
    if (read_narr_parameter_values (layers, "HGT_1", hgt1) != SUCCESS)
    {
        RETURN_ERROR ("Failed loading HGT_1 parameters", FUNC_NAME, FAILURE);
    }

    /* Read in NARR specific humidity for time before Landsat acqusition */
    if (read_narr_parameter_values (layers, "SPFH_1", spfh1) != SUCCESS)
    {
        RETURN_ERROR ("Failed loading SPFH_1 parameters", FUNC_NAME, FAILURE);
    }

    /* Read in NARR temperature for time before Landsat acqusition */
    if (read_narr_parameter_values (layers, "TMP_1", tmp1) != SUCCESS)
    {
        RETURN_ERROR ("Failed loading TMP_1 parameters", FUNC_NAME, FAILURE);
    }

    /* Read in NARR height for time after Landsat acqusition */
    if (read_narr_parameter_values (layers, "HGT_2", hgt2) != SUCCESS)
    {
        RETURN_ERROR ("Failed loading HGT_2 parameters", FUNC_NAME, FAILURE);
    }

    /* Read in NARR specific humidity for time after Landsat acqusition */
    if (read_narr_parameter_values (layers, "SPFH_2", spfh2) != SUCCESS)
    {
        RETURN_ERROR ("Failed loading SPFH_2 parameters", FUNC_NAME, FAILURE);
    }

    /* Read in NARR temperature for time after Landsat acqusition */
    if (read_narr_parameter_values (layers, "TMP_2", tmp2) != SUCCESS)
    {
        RETURN_ERROR ("Failed loading TMP_2 parameters", FUNC_NAME, FAILURE);
    }

    /* ==================================================================== */

    narr_hgt1 = (double **) allocate_2d_array (P_LAYER, num_points,
                                              sizeof (double));
    if (narr_hgt1 == NULL)
    {
        RETURN_ERROR ("Allocating narr_hgt1 memory", FUNC_NAME, FAILURE);
    }

    narr_hgt2 = (double **) allocate_2d_array (P_LAYER, num_points,
                                              sizeof (double));
    if (narr_hgt2 == NULL)
    {
        RETURN_ERROR ("Allocating narr_hgt2 memory", FUNC_NAME, FAILURE);
    }

    /* Allocate memory for humidity of NARR points within the rectangular */
    narr_spfh1 = (double **) allocate_2d_array (P_LAYER, num_points,
                                               sizeof (double));
    if (narr_spfh1 == NULL)
    {
        RETURN_ERROR ("Allocating narr_spfh1 memory", FUNC_NAME, FAILURE);
    }

    narr_spfh2 = (double **) allocate_2d_array (P_LAYER, num_points,
                                               sizeof (double));
    if (narr_spfh2 == NULL)
    {
        RETURN_ERROR ("Allocating narr_spfh2 memory", FUNC_NAME, FAILURE);
    }

    /* Allocate memory for temperature of NARR points within the rectangular */
    narr_tmp1 = (double **) allocate_2d_array (P_LAYER, num_points,
                                              sizeof (double));
    if (narr_tmp1 == NULL)
    {
        RETURN_ERROR ("Allocating narr_tmp1 memory", FUNC_NAME, FAILURE);
    }

    narr_tmp2 = (double **) allocate_2d_array (P_LAYER, num_points,
                                              sizeof (double));
    if (narr_tmp2 == NULL)
    {
        RETURN_ERROR ("Allocating narr_tmp2 memory", FUNC_NAME, FAILURE);
    }

    /* ==================================================================== */

    int index;
    for (layer = 0; layer < P_LAYER; layer++)
    {
        for (row = min_row; row <= max_row; row++)
        {
            for (col = min_col; col <= max_col; col++)
            {
                index = (row - min_row) * num_cols + (col - min_col);

                narr_hgt1[layer][index] = hgt1[layer][row * NARR_COLS + col];
                narr_spfh1[layer][index] = spfh1[layer][row * NARR_COLS + col];
                narr_tmp1[layer][index] = tmp1[layer][row * NARR_COLS + col];
                narr_hgt2[layer][index] = hgt2[layer][row * NARR_COLS + col];
                narr_spfh2[layer][index] = spfh2[layer][row * NARR_COLS + col];
                narr_tmp2[layer][index] = tmp2[layer][row * NARR_COLS + col];
            }
        }
    }

    /* ==================================================================== */

    /* Release memory */
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
    pressure = (double **) allocate_2d_array (P_LAYER, num_points,
                                             sizeof (double));
    if (pressure == NULL)
    {
        RETURN_ERROR ("Allocating pressure memory", FUNC_NAME, FAILURE);
    }

    for (layer = 0; layer < P_LAYER; layer++)
    {
        for (point = 0; point < num_points; point++)
        {
            pressure[layer][point] = layers[layer];
        }
    }

    /* ==================================================================== */

    narr_height1 = (double **) allocate_2d_array (P_LAYER, num_points,
                                                 sizeof (double));
    if (narr_height1 == NULL)
    {
        RETURN_ERROR ("Allocating narr_height1 memory", FUNC_NAME, FAILURE);
    }

    narr_height2 = (double **) allocate_2d_array (P_LAYER, num_points,
                                                 sizeof (double));
    if (narr_height2 == NULL)
    {
        RETURN_ERROR ("Allocating narr_height2 memory", FUNC_NAME, FAILURE);
    }

    narr_rh1 = (double **) allocate_2d_array (P_LAYER, num_points,
                                             sizeof (double));
    if (narr_rh1 == NULL)
    {
        RETURN_ERROR ("Allocating narr_rh1 memory", FUNC_NAME, FAILURE);
    }

    narr_rh2 = (double **) allocate_2d_array (P_LAYER, num_points,
                                             sizeof (double));
    if (narr_rh2 == NULL)
    {
        RETURN_ERROR ("Allocating narr_rh2 memory", FUNC_NAME, FAILURE);
    }

    /* ==================================================================== */

    /* convert grib data to variables to be input to MODTRAN */
    if (convert_geopotential_geometric (num_points, points->lat,
                                        narr_hgt1, narr_height1) != SUCCESS)
    {
        RETURN_ERROR ("Calling convert_geopotential_geometric for height 1",
                      FUNC_NAME, FAILURE);
    }

    if (convert_geopotential_geometric (num_points, points->lat,
                                        narr_hgt2, narr_height2) != SUCCESS)
    {
        RETURN_ERROR ("Calling convert_geopotential_geometric for height 2",
                      FUNC_NAME, FAILURE);
    }

    if (convert_sh_rh (num_points, points->lat, narr_spfh1, narr_tmp1,
                       pressure, narr_rh1) != SUCCESS)
    {
        RETURN_ERROR ("Calling convert_sh_rh for temp 1", FUNC_NAME,
                      FAILURE);
    }

    if (convert_sh_rh (num_points, points->lat, narr_spfh2, narr_tmp2,
                       pressure, narr_rh2) != SUCCESS)
    {
        RETURN_ERROR ("Calling convert_sh_rh for temp 2", FUNC_NAME,
                      FAILURE);
    }

    /* ==================================================================== */

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
    hour1 = (double) (input->meta.acq_date.hour - rem1);
    hour2 = (double) (input->meta.acq_date.hour + rem2);
    inv_hour_diff = 1.0 / (hour2 - hour1);

    /* Round to the nearest minute */
    if ((input->meta.acq_date.second - 30.0) >= MINSIGMA)
        input->meta.acq_date.minute++;

    /* convert hour-min acquisition time to decimal time */
    time = (double) input->meta.acq_date.hour
           + (double) input->meta.acq_date.minute / 60.0F;
    time_diff = time - hour1;

    /* ==================================================================== */

    /* Allocate memory */
    narr_height = (double **) allocate_2d_array (P_LAYER, num_points,
                                                sizeof (double));
    if (narr_height == NULL)
    {
        RETURN_ERROR ("Allocating narr_height memory", FUNC_NAME, FAILURE);
    }

    /* Allocate memory */
    narr_rh = (double **) allocate_2d_array (P_LAYER, num_points,
                                            sizeof (double));
    if (narr_rh == NULL)
    {
        RETURN_ERROR ("Allocating narr_rh memory", FUNC_NAME, FAILURE);
    }

    /* Allocate memory */
    narr_tmp = (double **) allocate_2d_array (P_LAYER, num_points,
                                             sizeof (double));
    if (narr_tmp == NULL)
    {
        RETURN_ERROR ("Allocating narr_tmp memory", FUNC_NAME, FAILURE);
    }

    /* linearly interpolate geometric height, relative humidity, and
       temperature for NARR points within Landsat scene this is the NARR data
       corresponding to the acquisition time of the Landsat image converted to
       appropriated variable for MODTRAN input */
    for (layer = 0; layer < P_LAYER; layer++)
    {
        for (point = 0; point < num_points; point++)
        {
            narr_height[layer][point] =
                narr_height1[layer][point] + time_diff
                * ((narr_height2[layer][point] - narr_height1[layer][point])
                   * inv_hour_diff);
            narr_rh[layer][point] =
                narr_rh1[layer][point] + time_diff
                * ((narr_rh2[layer][point] - narr_rh1[layer][point])
                   * inv_hour_diff);
            narr_tmp[layer][point] =
                narr_tmp1[layer][point] + time_diff
                * ((narr_tmp2[layer][point] - narr_tmp1[layer][point])
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
    stan_height = (double *) malloc (STANDARD_LAYERS * sizeof (double));
    if (stan_height == NULL)
    {
        RETURN_ERROR ("Allocating stan_height memory", FUNC_NAME, FAILURE);
    }

    stan_pre = (double *) malloc (STANDARD_LAYERS * sizeof (double));
    if (stan_pre == NULL)
    {
        RETURN_ERROR ("Allocating stan_pre memory", FUNC_NAME, FAILURE);
    }

    stan_temp = (double *) malloc (STANDARD_LAYERS * sizeof (double));
    if (stan_temp == NULL)
    {
        RETURN_ERROR ("Allocating stan_temp memory", FUNC_NAME, FAILURE);
    }

    stan_rh = (double *) malloc (STANDARD_LAYERS * sizeof (double));
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
    points->num_modtran_runs = num_modtran_runs;

    /* Allocate memory */
    points->modtran_runs =
        (MODTRAN_INFO *) malloc (num_modtran_runs * sizeof (MODTRAN_INFO));
    if (points->modtran_runs == NULL)
    {
        RETURN_ERROR ("Allocating modtran_runs memory", FUNC_NAME, FAILURE);
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

    char *dummy = NULL;
    dummy = getcwd (curr_path, PATH_MAX);
    if (dummy == NULL)
    {
        RETURN_ERROR ("Retrieving current working directory",
                      FUNC_NAME, FAILURE);
    }

    /* Allocate some temp memory */
    temp_height = (double *) malloc (MAX_MODTRAN_LAYER * sizeof (double));
    if (temp_height == NULL)
    {
        RETURN_ERROR ("Allocating temp_height memory", FUNC_NAME, FAILURE);
    }

    temp_pressure = (double *) malloc (MAX_MODTRAN_LAYER * sizeof (double));
    if (temp_pressure == NULL)
    {
        RETURN_ERROR ("Allocating temp_pressure memory", FUNC_NAME, FAILURE);
    }

    temp_temp = (double *) malloc (MAX_MODTRAN_LAYER * sizeof (double));
    if (temp_temp == NULL)
    {
        RETURN_ERROR ("Allocating temp_temp memory", FUNC_NAME, FAILURE);
    }

    temp_rh = (double *) malloc (MAX_MODTRAN_LAYER * sizeof (double));
    if (temp_rh == NULL)
    {
        RETURN_ERROR ("Allocating temp_rh memory", FUNC_NAME, FAILURE);
    }

    /* Create a point list / directory names that can be used later
       to delete them */
    point_list_fd = fopen ("point_list.txt", "w");
    if (point_list_fd == NULL)
    {
        RETURN_ERROR ("Opening file: point_list.txt\n", FUNC_NAME, FAILURE);
    }

    for (point = 0; point < num_points; point++)
    {
        /* ****************************************************************
           MODTRAN uses longitudinal degree values from 0 to 360 starting at
           Greenwich and moving west.

           The following logic fixes the longitude to be for MODTRAN and we
           also use those values for generating the output directory names.
           **************************************************************** */
        if (points->lon[point] < 0)
        {
            /* We are a west longitude value so negate it to the positive
               equivalent value.
               i.e.   W40 normally represented as -40 would be changed to a
                      positive 40 value. */
            points->lon[point] = -points->lon[point];
        }
        else
        {
            /* We are a east longitude value so fix it to be greater than
               180 west.
               i.e.   E10 would be turned into W350, and be positive 350 not
                      negative.  */
            points->lon[point] = 360.0 - points->lon[point];
        }

        /* Figure out the lat and lon strings to use.
           MODTRAN tape files are finicky about value locations and size,
           so the following adjusts for the values less than 100. */
        if (points->lat[point] < 100.0)
            snprintf (lat_str, sizeof (lat_str), "%6.3f", points->lat[point]);
        else
            snprintf (lat_str, sizeof (lat_str), "%6.2f", points->lat[point]);

        if (points->lon[point] < 100.0)
            snprintf (lon_str, sizeof (lon_str), "%6.3f", points->lon[point]);
        else
            snprintf (lon_str, sizeof (lon_str), "%6.2f", points->lon[point]);

        /* Create the name of the directory for the current NARR point */
        snprintf (current_point, sizeof (current_point),
                  "%s_%s", lat_str, lon_str);

        /* Create the directory */
#if 0
        snprintf(msg_str, sizeof (msg_str),
                 "Creating directory [%s]", current_point);
        LOG_MESSAGE (msg_str, FUNC_NAME);
#endif
        if (mkdir (current_point, 0755) != SUCCESS)
        {
            if (errno != EEXIST)
            {
                snprintf (msg_str, sizeof (msg_str),
                          "Failed creating directory [%s]", current_point);
                RETURN_ERROR (msg_str, FUNC_NAME, FAILURE);
            }
        }

        /* Save the directory name to the point list file */
        fprintf (point_list_fd, "%s\n", current_point);

        /* determine latitude and longitude of current NARR point and insert
           into tail file */
        /* TODO - The sed commands could be combined into one to save
                  traversing the file three times.  But the file is small. */
        snprintf (command, sizeof (command),
                  "cat %s/modtran_tail.txt"
                  " | sed 's/latitu/%s/'"
                  " | sed 's/longit/%s/'"
                  " | sed 's/jay/%d/' > new_tail.txt",
                  lst_data_dir, lat_str, lon_str, input->meta.acq_date.doy);
        if (system (command) != SUCCESS)
        {
            RETURN_ERROR ("Failed creating new_tail.txt", FUNC_NAME, FAILURE);
        }

        /* Clear the temp memory */
        for (layer = 0; layer < MAX_MODTRAN_LAYER; layer++)
        {
            temp_height[layer] = 0.0;
            temp_pressure[layer] = 0.0;
            temp_temp[layer] = 0.0;
            temp_rh[layer] = 0.0;
        }

        /* set lowest altitude is the first geometric height at that NARR
           point (if positive) and (if negative set to zero) */
        if (narr_height[0][point] < 0)
            gndalt[0] = 0.0;
        else
            gndalt[0] = narr_height[0][point];

        /* iterate through all ground altitudes at which MODTRAN is run */
        for (elevation = 0; elevation < NUM_ELEVATIONS; elevation++)
        {
            /* create a directory for the current height */
            snprintf (current_gdalt, sizeof (current_gdalt),
                      "%s/%5.3f", current_point, gndalt[elevation]);

            /* Create the directory */
#if 0
            snprintf (msg_str, sizeof (msg_str),
                      "Creating directory [%s]", current_gdalt);
            LOG_MESSAGE (msg_str, FUNC_NAME);
#endif
            if (mkdir (current_gdalt, 0755) != SUCCESS)
            {
                if (errno != EEXIST)
                {
                    snprintf (msg_str, sizeof (msg_str),
                              "Failed creating directory [%s]", current_gdalt);
                    RETURN_ERROR (msg_str, FUNC_NAME, FAILURE);
                }
            }

            /* determine layers below current gndalt and closest layer
               above and below */
            for (layer = 0; layer < P_LAYER; layer++)
            {
                if ((narr_height[layer][point] - gndalt[elevation]) >= MINSIGMA)
                {
                    layer_below = layer - 1;
                    layer_above = layer;
                    break;
                }
            }

            if (layer_below < 0)
            {
                layer_below = 0;
                layer_above = 1;
            }

            /* To save divisions */
            inv_height_diff = 1.0 / (narr_height[layer_above][point]
                                     - narr_height[layer_below][point]);

            /* linearly interpolate pressure, temperature, and relative 
               humidity to gndalt for lowest layer */
            new_pressure = pressure[layer_below][point]
                           + (gndalt[elevation]
                              - narr_height[layer_below][point])
                           * ((pressure[layer_above][point]
                               - pressure[layer_below][point])
                              * inv_height_diff);
            new_temp = narr_tmp[layer_below][point]
                       + (gndalt[elevation] - narr_height[layer_below][point])
                       * ((narr_tmp[layer_above][point]
                           - narr_tmp[layer_below][point])
                          * inv_height_diff);
            new_rh = narr_rh[layer_below][point]
                     + (gndalt[elevation] - narr_height[layer_below][point])
                     * ((narr_rh[layer_above][point]
                         - narr_rh[layer_below][point])
                        * inv_height_diff);

            /* create arrays containing only layers to be included in current
               tape5 file */
            curr_layer = 0;
            temp_height[curr_layer] = gndalt[elevation];
            temp_pressure[curr_layer] = new_pressure;
            temp_temp[curr_layer] = new_temp;
            temp_rh[curr_layer] = new_rh;
            curr_layer++;

            for (layer = layer_above; layer < P_LAYER; layer++)
            {
                temp_height[curr_layer] = narr_height[layer][point];
                temp_pressure[curr_layer] = pressure[layer][point];
                temp_temp[curr_layer] = narr_tmp[layer][point];
                temp_rh[curr_layer] = narr_rh[layer][point];
                curr_layer++;
            }


            /* MODTRAN throws an error when there are two identical layers in
               the tape5 file, if the current ground altitude and the next
               highest layer are close enough, eliminate interpolated layer */
            if (fabs (gndalt[elevation] - narr_height[layer_above][point])
                < 0.001)
            {
                curr_layer = 0;
                for (layer = layer_above; layer < P_LAYER; layer++)
                {
                    temp_height[curr_layer] = narr_height[layer][point];
                    temp_pressure[curr_layer] = pressure[layer][point];
                    temp_temp[curr_layer] = narr_tmp[layer][point];
                    temp_rh[curr_layer] = narr_rh[layer][point];
                    curr_layer++;
                }
            }

            /* determine maximum height of NARR layers and where the standard 
               atmosphere is greater than this */
            std_layer = 0;
            for (layer = 0; layer < STANDARD_LAYERS; layer++)
            {
                if (stan_height[layer] > narr_height[P_LAYER - 1][point])
                {
                    counter[std_layer] = layer;
                    std_layer++;
                }
            }

            /* If there are more than 2 layers above the highest NARR layer,
               then we need to interpolate a value between the highest NARR
               layer and the 2nd standard atmosphere layer above the NARR
               layers to create a smooth transition between the NARR layers
               and the standard upper atmosphere */
            if (std_layer >= 3)
            {
#if 0
                snprintf (msg_str, sizeof (msg_str),
                          "Adding interpolated layer between the NARR layers"
                          " and the standard atmosphere.");
                LOG_MESSAGE (msg_str, FUNC_NAME);
#endif

                /* To save divisions */
                inv_height_diff = 1.0 / (stan_height[counter[2]]
                                         - temp_height[curr_layer - 1]);

                new_height = (stan_height[counter[2]]
                              + temp_height[curr_layer - 1]) / 2.0;
                new_pressure = temp_pressure[curr_layer - 1]
                               + (new_height - temp_height[curr_layer - 1])
                               * ((stan_pre[counter[2]]
                                   - temp_pressure[curr_layer - 1])
                                  * inv_height_diff);
                new_temp = temp_temp[curr_layer - 1]
                           + (new_height - temp_height[curr_layer - 1])
                           * ((stan_temp[counter[2]]
                               - temp_temp[curr_layer - 1])
                              * inv_height_diff);
                new_rh = temp_rh[curr_layer - 1]
                         + (new_height - temp_height[curr_layer - 1])
                         * ((stan_rh[counter[2]] - temp_rh[curr_layer - 1])
                            * inv_height_diff);

                /* concatenate NARR layers, new layer, and standard atmosphere
                   layers */
                temp_height[curr_layer] = new_height;
                temp_pressure[curr_layer] = new_pressure;
                temp_temp[curr_layer] = new_temp;
                temp_rh[curr_layer] = new_rh;
                curr_layer++;
            }

            /* Add the remaining standard atmosphere layers */
            for (layer = 2; layer < std_layer; layer++)
            {
                temp_height[curr_layer] = stan_height[counter[layer]];
                temp_pressure[curr_layer] = stan_pre[counter[layer]];
                temp_temp[curr_layer] = stan_temp[counter[layer]];
                temp_rh[curr_layer] = stan_rh[counter[layer]];
                curr_layer++;
            }

            /* write atmospheric layers to a text file in format proper for
               tape5 file */
            fd = fopen ("temp_layers.txt", "w");
            if (fd == NULL)
            {
                RETURN_ERROR ("Opening file: temp_layers.txt\n",
                              FUNC_NAME, FAILURE);
            }

            /* Write out the intermediate file */
            for (layer = 0; layer < curr_layer; layer++)
            {
                fprintf (fd, "%10.3f%10.3e%10.3e%10.3e%10.3e%10.3e%16s\n",
                         temp_height[layer], temp_pressure[layer],
                         temp_temp[layer], temp_rh[layer], 0.0, 0.0,
                         "AAH             ");
            }

            /* Close the intermediate file */
            if (fclose (fd) != SUCCESS)
            {
                RETURN_ERROR ("Closing file: temp_layers.txt\n",
                              FUNC_NAME, FAILURE);
            }
            fd = NULL;

            /* determine number of layers for current ground altitude and
               insert into head file */
            snprintf (command, sizeof (command),
                      "cat %s/modtran_head.txt"
                      " | sed 's/nml/%d/'"
                      " | sed 's/gdalt/%5.3f/'"
                      " > base_head.txt",
                      lst_data_dir, curr_layer, gndalt[elevation]);
            if (system (command) != SUCCESS)
            {
                RETURN_ERROR ("Failed creating tempHead.txt", FUNC_NAME,
                              FAILURE);
            }

            /* iterate through [temperature,albedo] pairs at which to run
               MODTRAN */
            for (temperature = 0; temperature <= 2; temperature++)
            {
                /* create directory for the current temperature */
                snprintf (current_temp, sizeof (current_temp),
                          "%s/%s", current_gdalt, temp_strs[temperature]);
#if 0
                snprintf (msg_str, sizeof (msg_str),
                          "Creating directory [%s]", current_temp);
                LOG_MESSAGE (msg_str, FUNC_NAME);
#endif
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

                /* create directory for the current albedo */
                snprintf (current_alb, sizeof (current_alb),
                          "%s/%3.1f", current_temp, alb[temperature]);
#if 0
                snprintf (msg_str, sizeof (msg_str),
                          "Creating directory [%s]", current_alb);
                LOG_MESSAGE (msg_str, FUNC_NAME);
#endif
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

                /* Concatenate the base head and atmospheric layers and tail
                   files to create a tape 5 file for MODTRAN specific to this
                   location and ground altitude with variables for temperature
                   and albedo substituted file */
                snprintf (command, sizeof (command),
                          "cat base_head.txt"
                          " temp_layers.txt"
                          " new_tail.txt"
                          " | sed 's/tmp/%s/'"
                          " | sed 's/alb/%4.2f/'"
                          " > %s/tape5",
                          temp_strs[temperature], alb[temperature],
                          current_alb);
                if (system (command) != SUCCESS)
                {
                    RETURN_ERROR ("Failed creating tape5", FUNC_NAME, FAILURE);
                }

                /* create string for case list containing the location of the
                   current tape5 file

                   create string for command list containing the commands for
                   the MODTRAN run

                   iterate entry count */
                case_counter = point * NUM_ELEVATIONS * 3
                               + elevation * 3 + temperature;
                snprintf (points->modtran_runs[case_counter].path, PATH_MAX,
                          "%s", current_alb);
                snprintf (points->modtran_runs[case_counter].command, PATH_MAX,
                          "cd %s; ln -s %s DATA; %s/modtran",
                          points->modtran_runs[case_counter].path,
                          modtran_data_dir, modtran_path);

                points->modtran_runs[case_counter].latitude =
                    points->lat[point];
                points->modtran_runs[case_counter].longitude =
                    points->lon[point];
                points->modtran_runs[case_counter].height =
                    gndalt[elevation];
            } /* END - Temperature Albedo Pairs */
        } /* END - ground altitude ran by MODTRAN */
    } /* END - number of points */

    /* Close the point_list.txt file */
    if (fclose (point_list_fd) != SUCCESS)
    {
        RETURN_ERROR ("Closing file: point_list.txt\n", FUNC_NAME, FAILURE);
    }

    /* Free the temp memory */
    free(temp_height);
    free(temp_pressure);
    free(temp_temp);
    free(temp_rh);
    temp_height = NULL;
    temp_pressure = NULL;
    temp_temp = NULL;
    temp_rh = NULL;

    /* Free the standard atmosphere memory */
    free(stan_height);
    free(stan_pre);
    free(stan_temp);
    free(stan_rh);
    stan_height = NULL;
    stan_pre = NULL;
    stan_temp = NULL;
    stan_rh = NULL;

    if (debug)
    {
        /* write case_list.txt to a file */
        fd = fopen ("case_list.txt", "w");
        if (fd == NULL)
        {
            RETURN_ERROR ("Opening file: case_list.txt\n",
                          FUNC_NAME, FAILURE);
        }

        /* Write out the case_list.txt file */
        for (index = 0; index < num_points * NUM_ELEVATIONS * 3; index++)
        {
            fprintf (fd, "%s\n", points->modtran_runs[index].path);
        }

        /* Close the case_list.txt file */
        if (fclose (fd) != SUCCESS)
        {
            RETURN_ERROR ("Closing file: case_list.txt\n",
                          FUNC_NAME, FAILURE);
        }

        /* write command_list.txt to a file */
        fd = fopen ("command_list.txt", "w");
        if (fd == NULL)
        {
            RETURN_ERROR ("Opening file: command_list.txt\n",
                          FUNC_NAME, FAILURE);
        }

        /* Write out the command_list.txt file */
        for (index = 0; index < num_modtran_runs; index++)
        {
            fprintf (fd, "%s\n", points->modtran_runs[index].command);
        }

        /* Close the command_list.txt file */
        if (fclose (fd) != SUCCESS)
        {
            RETURN_ERROR ("Closing file: command_list.txt\n",
                          FUNC_NAME, FAILURE);
        }
        fd = NULL;
    }

    /* Free remaining memory allocations */
    if (free_2d_array ((void **) pressure) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: pressure\n", FUNC_NAME, FAILURE);
    }
    pressure = NULL;

    if (free_2d_array ((void **) narr_height) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: narr_height\n", FUNC_NAME, FAILURE);
    }
    narr_height = NULL;

    if (free_2d_array ((void **) narr_rh) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: narr_rh\n", FUNC_NAME, FAILURE);
    }
    narr_rh = NULL;

    if (free_2d_array ((void **) narr_tmp) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: narr_tmp\n", FUNC_NAME, FAILURE);
    }
    narr_tmp = NULL;

    return SUCCESS;
}
