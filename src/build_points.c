
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
    double **lat,
    double **lon
)
{
    char FUNC_NAME[] = "read_narr_coordinates";
    int row;
    int col;
    int count;
    int grid_row;
    int grid_col;
    double grid_lat;
    double grid_lon;
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
    for (row = 0; row < NARR_ROWS; row++)
    {
        for (col = 0; col < NARR_COLS; col++)
        {
            /* File Format:
               Grid_Column Grid_Row Grid_Latitude Grid_Longitude
             */
            if (fscanf (fd, "%d %d %lf %lf",
                        &grid_col, &grid_row, &grid_lat, &grid_lon)
                == EOF)
            {
                RETURN_ERROR ("End of file (EOF) is met before"
                              " NARR_ROWS * NARR_COLS lines",
                              FUNC_NAME, FAILURE);
            }

            lat[row][col] = grid_lat;

            /* TODO - Should think about fixing the input file, so that this
                      confusing conversion is not needed.

                      When you read the file data, it is as if you are
                      reading the values from the lower left to the upper
                      right as applied to the earth.  And the values being
                      read in start with an origin somewhere around the lower
                      left, hence the need for the following conversion.

               NOTE - If this is changed here, then else-where in the code
                      will break. */
            if (grid_lon > 180.0)
                lon[row][col] = 360.0 - grid_lon;
            else
                lon[row][col] = -grid_lon;
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
    Sep 2014    Song Guo         Original Development
    Feb 2015    Ron Dilley       Renamed and modified to use the
                                 REANALYSIS_POINTS data structure.
*****************************************************************************/
void convert_ll_to_utm
(
    Input_t *input,
    REANALYSIS_POINTS *points /* I/O: The coordinate points to be used */
)
{
    int point;
    double a = UTM_EQUATORIAL_RADIUS; /* equatorial radius */
    double b = UTM_POLAR_RADIUS;      /* polar radius */
    double k0 = UTM_SCALE_FACTOR;     /* scale factor */
    double ecc;                       /* eccentricity */
    double ecc_prime_sqrd;            /* prime of eccentricity squared */
    double n;
    double nu;
    double a0;
    double b0;
    double c0;
    double d0;
    double e0;
    double ki;
    double kii;
    double kiii;
    double kiv;
    double kv;
    double zone_cm;
    double p;
    double lat_rad;
    double s;

    /* Calculate zone central meridian in degrees and radians */
    zone_cm = (6.0 * input->meta.zone) - 183.0;

    ecc = sqrt (1.0 - ((b / a) * (b / a)));
    ecc_prime_sqrd = (ecc * ecc) / (1.0 - ecc * ecc);
    n = (a - b) / (a + b);

    /* Calculate meridional arc length */
    a0 = a * (1.0 - n
              + (5.0 * n * n / 4.0) * (1.0 - n)
              + (81.0 * pow (n, 4) / 64.0) * (1.0 - n));

    b0 = (3.0 * a * n / 2.0) * (1.0 - n
                                - (7.0 * n * n / 8.0) * (1.0 - n)
                                + 55.0 * pow (n, 4) / 64.0);

    c0 = (15.0 * a * n * n / 16.0) * (1.0 - n
                                      + (3.0 * n * n / 4.0) * (1.0 - n));

    d0 = (35 * a * pow (n, 3) / 48.0) * (1.0 - n
                                         + 11.0 * n * n / 16.0);

    e0 = (315.0 * a * pow (n, 4) / 51.0) * (1.0 - n);

    for (point = 0; point < points->num_points; point++)
    {
        /* Distance from the point to the zone central meridian in radians */
        p = (points->lon[point] - zone_cm) * RADIANS_PER_DEGREE;

        /* Convert latitude for the point to radians */
        lat_rad = points->lat[point] * RADIANS_PER_DEGREE;

        nu = a / sqrt (1.0 - (ecc * sin (lat_rad))
                             * (ecc * sin (lat_rad)));

        s = a0 * lat_rad
            - b0 * sin (2 * lat_rad)
            + c0 * sin (4 * lat_rad)
            - d0 * sin (6 * lat_rad)
            + e0 * sin (8 * lat_rad);

        /* Coefficients for UTM coordinates */
        ki = s * k0;

        kii = nu * sin (lat_rad) * cos (lat_rad) * k0 / 2.0;

        kiii = (nu * sin (lat_rad) * pow (cos (lat_rad), 3) / 24.0)
               * (5.0 - tan (lat_rad) * tan (lat_rad)
                  + 9.0 * ecc_prime_sqrd * pow (cos (lat_rad), 2)
                  + 4.0 * ecc_prime_sqrd * ecc_prime_sqrd
                    * pow (cos (lat_rad), 4))
               * k0;

        kiv = nu * cos (lat_rad) * k0;

        kv = pow (cos (lat_rad), 3)
             * (nu / 6.0)
             * (1.0 - tan (lat_rad) * tan (lat_rad)
                + ecc_prime_sqrd * cos (lat_rad) * cos (lat_rad))
             * k0;

        /* Calculate UTM coordinates */
        points->utm_easting[point] = UTM_FALSE_EASTING
                                     + (kiv * p + kv * pow (p, 3));
        points->utm_northing[point] = (ki + kii * p * p + kiii * pow (p, 4));
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
    char msg[PATH_MAX];

    double **lat;
    double **lon;

    int row;
    int col;
    int min_row;
    int max_row;
    int min_col;
    int max_col;

    int num_bytes;
    int index;

    double buffered_ul_lat;
    double buffered_ul_lon;
    double buffered_lr_lat;
    double buffered_lr_lon;

    /* Grab the environment path to the LST_DATA_DIR */
    lst_data_dir = getenv ("LST_DATA_DIR");
    if (lst_data_dir == NULL)
    {
        RETURN_ERROR ("LST_DATA_DIR environment variable is not set",
                      FUNC_NAME, FAILURE);
    }

    /* Allocate 2d memory to hold the coordinates */
    lat = (double **) allocate_2d_array (NARR_ROWS, NARR_COLS, sizeof (double));
    if (lat == NULL)
    {
        RETURN_ERROR ("Allocating lat memory", FUNC_NAME, FAILURE);
    }

    lon = (double **) allocate_2d_array (NARR_ROWS, NARR_COLS, sizeof (double));
    if (lon == NULL)
    {
        RETURN_ERROR ("Allocating lon memory", FUNC_NAME, FAILURE);
    }

    /* Read the coordinates into memory */
    if (read_narr_coordinates (lst_data_dir, lat, lon) != SUCCESS)
    {
        RETURN_ERROR ("Failed loading HGT_1 parameters", FUNC_NAME, FAILURE);
    }

    /* expand range to include NARR points outside image for edge pixels */
    /* TODO - 0.2deg at higher latitudes will not be sufficient for the
              longitudinal test, since at lat(72deg) 1deg lon = 34504.22meters
              and the NARR data is 34k between points.

              This is probably only a CONUS quick and dirty solution.

       NOTE - MERRA is even farther apart so this will not work for that. */
    buffered_ul_lat = input->meta.ul_geo_corner.lat + 0.2;
    buffered_ul_lon = input->meta.ul_geo_corner.lon - 0.2;
    buffered_lr_lat = input->meta.lr_geo_corner.lat - 0.2;
    buffered_lr_lon = input->meta.lr_geo_corner.lon + 0.2;

    /* determine what points in the NARR dataset fall within our buffered
       Landsat area using logical operators lessThanLat and greaterThanLat
       are values where the NARR values are less than or greater than the
       edges of the Landsat corners values respectively pixels that are true
       in both fall within the Landsat scene the same thing is done with
       longitude values */
    min_row = 1000;
    max_row = 0;
    min_col = 1000;
    max_col = 0;
    for (row = 0; row < NARR_ROWS; row++)
    {
        for (col = 0; col < NARR_COLS; col++)
        {
            if ((buffered_ul_lat > lat[row][col])
                && (buffered_lr_lat < lat[row][col])
                && (buffered_ul_lon < lon[row][col])
                && (buffered_lr_lon > lon[row][col]))
            {
                min_row = min (min_row, row);
                max_row = max (max_row, row);
                min_col = min (min_col, col);
                max_col = max (max_col, col);
            }
        }
    }

    /* Save these in the points structure */
    points->ul_lat = buffered_ul_lat;
    points->ul_lon = buffered_ul_lon;
    points->lr_lat = buffered_lr_lat;
    points->lr_lon = buffered_lr_lon;
    points->min_row = min_row;
    points->max_row = max_row;
    points->min_col = min_col;
    points->max_col = max_col;
    points->num_rows = max_row - min_row + 1;
    points->num_cols = max_col - min_col + 1;
    points->num_points = points->num_rows * points->num_cols;

    snprintf (msg, sizeof (msg), "min_row = %d\n", points->min_row);
    LOG_MESSAGE (msg, FUNC_NAME);

    snprintf (msg, sizeof (msg), "max_row = %d\n", points->max_row);
    LOG_MESSAGE (msg, FUNC_NAME);

    snprintf (msg, sizeof (msg), "min_col = %d\n", points->min_col);
    LOG_MESSAGE (msg, FUNC_NAME);

    snprintf (msg, sizeof (msg), "max_col = %d\n", points->max_col);
    LOG_MESSAGE (msg, FUNC_NAME);

    snprintf (msg, sizeof (msg), "num_rows = %d\n", points->num_rows);
    LOG_MESSAGE (msg, FUNC_NAME);

    snprintf (msg, sizeof (msg), "num_cols = %d\n", points->num_cols);
    LOG_MESSAGE (msg, FUNC_NAME);

    snprintf (msg, sizeof (msg), "num_points = %d\n", points->num_points);
    LOG_MESSAGE (msg, FUNC_NAME);

    /* Initialize this pointer */
    points->modtran_runs = NULL;

    /* Determine the number of bytes for the memory allocations */
    num_bytes = points->num_points * sizeof (double);

    /* Allocate memory for points within the rectangle */
    points->row = (double *) malloc (num_bytes);
    if (points->row == NULL)
    {
        RETURN_ERROR ("Allocating points row memory", FUNC_NAME, FAILURE);
    }

    points->col = (double *) malloc (num_bytes);
    if (points->col == NULL)
    {
        RETURN_ERROR ("Allocating points col memory", FUNC_NAME, FAILURE);
    }

    points->lat = (double *) malloc (num_bytes);
    if (points->lat == NULL)
    {
        RETURN_ERROR ("Allocating points lat memory", FUNC_NAME, FAILURE);
    }

    points->lon = (double *) malloc (num_bytes);
    if (points->lon == NULL)
    {
        RETURN_ERROR ("Allocating points lon memory", FUNC_NAME, FAILURE);
    }

    points->utm_easting = (double *) malloc (num_bytes);
    if (points->utm_easting == NULL)
    {
        RETURN_ERROR ("Allocating points utm_easting memory", FUNC_NAME,
                      FAILURE);
    }

    points->utm_northing = (double *) malloc (num_bytes);
    if (points->utm_northing == NULL)
    {
        RETURN_ERROR ("Allocating points utm_northing memory", FUNC_NAME,
                      FAILURE);
    }

    /* Retain only the points within the rectangle */
    for (row = min_row; row <= max_row; row++)
    {
        for (col = min_col; col <= max_col; col++)
        {
            index = (row - min_row) * points->num_cols + (col - min_col);

            /* Row and col are not used, but maintained for debugging */
            points->row[index] = row;
            points->col[index] = col;
            points->lat[index] = lat[row][col];
            points->lon[index] = lon[row][col];

            points->utm_easting[index] = 0.0;
            points->utm_northing[index] = 0.0;
        }
    }

    /* Convert lat/lon to UTM northing/easting*/
    convert_ll_to_utm (input, points);

    /* Free memory only used locally */
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
    free(points->row);
    free(points->col);
    free(points->lat);
    free(points->lon);
    free(points->utm_easting);
    free(points->utm_northing);

    points->modtran_runs = NULL;
    points->row = NULL;
    points->col = NULL;
    points->lat = NULL;
    points->lon = NULL;
    points->utm_easting = NULL;
    points->utm_northing = NULL;
}
