
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>


#include "const.h"
#include "utilities.h"
#include "2d_array.h"
#include "lst_types.h"
#include "input.h"


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
                      left, hence the need for the following conversion.
               NOTE - If this is changed here, the else-where in the code will
                      break. */
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


/*****************************************************************************
MODULE:  convert_ll_to_utm

PURPOSE: Convert latitude and longitude to UTM northing and easting

RETURN: NONE

HISTORY:
    Date        Programmer       Reason
    ----------  ---------------  ---------------------------------------------
    09/30/2014  Song Guo         Original Development
    Feb 2015    Ron Dilley       Renamed and modified to use the
                                 REANALYSIS_POINTS data structure.
*****************************************************************************/
void convert_ll_to_utm
(
    Input_t *input,
    REANALYSIS_POINTS *points /* I/O: The coordinate points to be used */
)
{
    int i;
    float a = 6378137.0;        /* equatorial radius */
    float b = 6356752.3;        /* polar radius */
    float k0 = 0.9996;          /* scale factor */
    float ecc;                  /* eccentricity */
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

    /* calculate zone central meridian in degrees and radians */
    zone_cm = (float) (6 * input->meta.zone - 183);
//    zone_cm_rad =  zone_cm * PI / 180.0;

    ecc = sqrt (1.0 - pow (b / a, 2));
    eprimesqrd = ecc * ecc / (1.0 - ecc * ecc);
    n = (a - b) / (a + b);

    /* calculate meridional arc length */
    a0 = a * (1.0 - n
              + (5.0 * n * n / 4.0) * (1.0 - n)
              + (81.0 * (float) pow (n, 4) / 64.0) * (1.0 - n));

    b0 = (3.0 * a * n / 2.0) * (1.0 - n
                                - (7.0 * n * n / 8.0) * (1.0 - n)
                                + 55.0 * (float) pow (n, 4) / 64.0);

    c0 = (15.0 * a * n * n / 16.0) * (1.0 - n
                                      + (3.0 * n * n / 4.0) * (1.0 - n));

    d0 = (35 * a * (float) pow (n, 3) / 48.0) * (1.0 - n
                                                 + 11.0 * n * n / 16.0);

    e0 = (315.0 * a * (float) pow (n, 4) / 51.0) * (1.0 - n);

    for (i = 0; i < points->num_points; i++)
    {
        delta_lon = points->lon[i] - zone_cm;
        p = delta_lon * PI / 180.0;

        /* convert lat and lon points from decimal degrees to radians */
        lat_rad = points->lat[i] * PI / 180.0;
//        lon_rad = lon[i]*PI/180.0;

//        rho = a*(1.0-e*e)/(pow((1.0-(e*sin(lat_rad))*(e*sin(lat_rad))),
//              (3.0/2.0)));
        nu = a / pow ((1.0 - (ecc * sin (lat_rad)) * (ecc * sin (lat_rad))),
                      (1.0 / 2.0));


        s = a0 * lat_rad
            - b0 * sin (2 * lat_rad)
            + c0 * sin (4 * lat_rad)
            - d0 * sin (6 * lat_rad)
            + e0 * sin (8 * lat_rad);

        /* coefficients for UTM coordinates */
        ki = s * k0;

        kii = nu * sin (lat_rad) * cos (lat_rad) * k0 / 2.0;

        kiii = (pow (nu * sin (lat_rad) * cos (lat_rad), 3) / 24.0)
               * pow ((5 - tan (lat_rad)), 2)
               + 9.0 * eprimesqrd * pow (cos (lat_rad), 2)
               + 4.0 * eprimesqrd * eprimesqrd * pow (cos (lat_rad), 4) * k0;

        kiv = nu * cos (lat_rad) * k0;

        kv = pow (cos (lat_rad), 3)
             * (nu / 6.0)
             * (1 - tan (lat_rad) * tan (lat_rad)
                + eprimesqrd * cos (lat_rad) * cos (lat_rad))
             * k0;

        /* calculate UTM coordinates */
        points->utm_easting[i] = 500000.0 + (kiv * p + kv * pow (p, 3));
        points->utm_northing[i] = (ki + kii * p * p + kiii * pow (p, 4));
    }
}


int build_points
(
    Input_t *input,           /* I: input structure */
    REANALYSIS_POINTS *points /* O: The coordinate points to be used */
)
{
    char FUNC_NAME[] = "build_points";
    char *lst_data_dir = NULL;

    int **eye;
    int **jay;
    float **lat;
    float **lon;

    int i;
    int j;
    int min_eye;
    int max_eye;
    int min_jay;
    int max_jay;

    int num_bytes;
    int index;

    float buffered_ul_lat;
    float buffered_ul_lon;
    float buffered_lr_lat;
    float buffered_lr_lon;

    /* Grab the environment path to the LST_DATA_DIR */
    lst_data_dir = getenv ("LST_DATA_DIR");
    if (lst_data_dir == NULL)
    {
        RETURN_ERROR ("LST_DATA_DIR environment variable is not set",
                      FUNC_NAME, FAILURE);
    }

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

    /* Read the coordinates into memory */
    if (read_narr_coordinates (lst_data_dir, eye, jay, lat, lon) != SUCCESS)
    {
        RETURN_ERROR ("Failed loading HGT_1 parameters", FUNC_NAME, FAILURE);
    }

    /* expand range to include NARR points outside image for edge pixels */
    /* TODO - 0.5deg at higher latitudes will not be sufficient for the
              longitudinal test, since at lat(72deg) 1deg lon = 34504.22meters
              and the NARR data is 34k between points.

              This is probably only a CONUS quick and dirty solution.

       NOTE - MERRA is even farther apart so this will not work for that. */
    buffered_ul_lat = input->meta.ul_geo_corner.lat + 0.5;
    buffered_ul_lon = input->meta.ul_geo_corner.lon - 0.5;
    buffered_lr_lat = input->meta.lr_geo_corner.lat - 0.5;
    buffered_lr_lon = input->meta.lr_geo_corner.lon + 0.5;

    /* determine what points in the NARR dataset fall within our buffered
       Landsat area using logical operators lessThanLat and greaterThanLat
       are values where the NARR values are less than or greater than the
       edges of the Landsat corners values respectively pixels that are true
       in both fall within the Landsat scene the same thing is done with
       longitude values */
    max_eye = 0;
    min_eye = 1000;
    max_jay = 0;
    min_jay = 1000;
    for (i = 0; i < NARR_ROW - 1; i++)
    {
        for (j = 0; j < NARR_COL - 1; j++)
        {
            if ((buffered_ul_lat > lat[i][j])
                && (buffered_lr_lat < lat[i][j])
                && (buffered_ul_lon < lon[i][j])
                && (buffered_lr_lon > lon[i][j]))
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

    /* Save these in the points structure */
    points->min_eye = min_eye;
    points->max_eye = max_eye;
    points->min_jay = min_jay;
    points->max_jay = max_jay;
    points->num_eyes = max_eye - min_eye + 1;
    points->num_jays = max_jay - min_jay + 1;
    points->num_points = points->num_eyes * points->num_jays;

    /* Determine the number of byte for the memory allocations */
    num_bytes = points->num_points * sizeof (float);

    /* Allocate memory for points within the rectangle */
    points->eye = (float *) malloc (num_bytes);
    if (points->eye == NULL)
    {
        RETURN_ERROR ("Allocating points eye memory", FUNC_NAME, FAILURE);
    }

    points->jay = (float *) malloc (num_bytes);
    if (points->jay == NULL)
    {
        RETURN_ERROR ("Allocating points jay memory", FUNC_NAME, FAILURE);
    }

    points->lat = (float *) malloc (num_bytes);
    if (points->lat == NULL)
    {
        RETURN_ERROR ("Allocating points lat memory", FUNC_NAME, FAILURE);
    }

    points->lon = (float *) malloc (num_bytes);
    if (points->lon == NULL)
    {
        RETURN_ERROR ("Allocating points lon memory", FUNC_NAME, FAILURE);
    }

    points->utm_easting = (float *) malloc (num_bytes);
    if (points->utm_easting == NULL)
    {
        RETURN_ERROR ("Allocating points utm_easting memory", FUNC_NAME,
                      FAILURE);
    }

    points->utm_northing = (float *) malloc (num_bytes);
    if (points->utm_northing == NULL)
    {
        RETURN_ERROR ("Allocating points utm_northing memory", FUNC_NAME,
                      FAILURE);
    }

    /* Retain only the points within the rectangle */
    for (j = min_jay; j <= max_jay; j++)
    {
        for (i = min_eye; i <= max_eye; i++)
        {
            index = (j - min_jay) * points->num_eyes + (i - min_eye);

            points->eye[index] = eye[i][j];
            points->jay[index] = jay[i][j];
            points->lat[index] = lat[i][j];
            points->lon[index] = lon[i][j];

            points->utm_easting[index] = 0.0;
            points->utm_northing[index] = 0.0;
        }
    }

    /* Convert lat/lon to UTM northing/easting*/
    convert_ll_to_utm (input, points);

    /* Free memory only used locally */
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

    return SUCCESS;
}


void free_points_memory
(
    REANALYSIS_POINTS *points /* I: The coordinate points */
)
{
    free(points->modtran_runs);
    free(points->eye);
    free(points->jay);
    free(points->lat);
    free(points->lon);
    free(points->utm_easting);
    free(points->utm_northing);

    points->modtran_runs = NULL;
    points->eye = NULL;
    points->jay = NULL;
    points->lat = NULL;
    points->lon = NULL;
    points->utm_easting = NULL;
    points->utm_northing = NULL;
}
