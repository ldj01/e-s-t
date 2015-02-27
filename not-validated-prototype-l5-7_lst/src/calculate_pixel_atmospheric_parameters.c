
#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <float.h>
#include <math.h>


#include "const.h"
#include "2d_array.h"
#include "utilities.h"
#include "input.h"
#include "output.h"
#include "lst_types.h"
#include "build_points.h"


#define NUM_PARAMETERS 3


/* Defines the distance to the current pixel, along with the index of the
   point
   So that we can find the index of the closest point to start determining the
   correct cell to use */
typedef struct
{
    int index;
    float distance;
} DISTANCE_ITEM;


/* Defines index locations in the vertices array for the current cell to be
   used for interpolation of the pixel */
typedef enum
{
    LL_POINT,
    UL_POINT,
    UR_POINT,
    LR_POINT,
    NUM_CELL_POINTS
} CELL_POINTS;


/* A qsort routine that can be used with the DISTANCE_ITEM items */
int qsort_distance_compare_function
(
    const void *distance_item_a,
    const void *distance_item_b
)
{
    float a = (*(DISTANCE_ITEM*)distance_item_a).distance;
    float b = (*(DISTANCE_ITEM*)distance_item_b).distance;

    if (a < b)
        return -1;
    else if (b < a)
        return 1;

    return 0;
}


/******************************************************************************
METHOD:  distance_in_utm

PURPOSE: Calculate distances between UTM coordiantes

RETURN: NONE

NOTE: SR(x) = (scale_factor / cos ((x - false_easting) / equatorial_radius))

NOTE: Simpson's Rule is applied for integrating the longitudinal distance
      from easting of first point to easting of second point.

      SR(x)dx ~= ((e2 - e0) / 6)
                 * (SR(e0) + 4 * SR((e0 + e2) / 2) + SR(e2))

      Where:
          e0 = easting of starting point
          e2 = easting of stopping point

******************************************************************************/
#define INV_UTM_EQUATORIAL_RADIUS (1.0 / UTM_EQUATORIAL_RADIUS)
#define INV_TWO (0.5) // (1.0 / 2.0)
#define INV_SIX (1.0 / 6.0)
void distance_in_utm
(
    float e0,
    float n0,
    float e2,
    float n2,
    float *distance
)
{
    /* The UTM coordinates we are using have the 500000 false easting applied
       to them, so we need to remove that before applying the distance
       calculation. */
    float e0_adj;
    float e1_term;
    float e2_adj;

    float sr_e0;
    float sr_e1;
    float sr_e2;

    float edist;

    e0_adj = e0 - UTM_FALSE_EASTING;
    e1_term = (e0 + e2) * INV_TWO;
    e2_adj = e2 - UTM_FALSE_EASTING;


    sr_e0 = UTM_SCALE_FACTOR / (cos (e0_adj * INV_UTM_EQUATORIAL_RADIUS));

    sr_e1 = UTM_SCALE_FACTOR / (cos (e1_term * INV_UTM_EQUATORIAL_RADIUS));

    sr_e2 = UTM_SCALE_FACTOR / (cos (e2_adj * INV_UTM_EQUATORIAL_RADIUS));

    edist = ((e2 - e0) * INV_SIX)
            * (sr_e0 + 4.0 * sr_e1 + sr_e2);

    *distance = sqrt (edist * edist + (n2 - n0) * (n2 - n0));
}


/******************************************************************************
METHOD:  interpolate_to_height

PURPOSE: Interpolate to height of current pixel

******************************************************************************/
void interpolate_to_height
(
    float **modtran_results, /* I: results from MODTRAN runs for a point */
    double interpolate_to,   /* I: current landsat pixel height */
    float *at_height         /* O: interpolated height for point */
)
{
    int i;
    int below = 0;
    int above = 0;

    float below_parameters[NUM_PARAMETERS];
    float above_parameters[NUM_PARAMETERS];

    float slope;
    float intercept;

    float above_height;
    float inv_height_diff; /* To remove the multiple divisions */

    /* Find the height to use that is below the interpolate_to height */
    for (i = 0; i < NUM_ELEVATIONS; i++)
    {
        if (modtran_results[i][LST_HEIGHT] < interpolate_to)
        {
            below = i; /* Last match will always be the one we want */
        }
    }

    /* Find the height to use that is equal to or above the interpolate_to
       height

       It will always be the same or the next height */
    above = below; /* Start with the same */
    if (above != (NUM_ELEVATIONS - 1))
    {
        /* Not the last height */

        /* Check to make sure that we are not less that the below height,
           indicating that our interpolate_to height is below the first
           height */
        if (! (interpolate_to < modtran_results[above][LST_HEIGHT]))
        {
            /* Use the next height, since it will be equal to or above our
               interpolate_to height */
            above++;
        }
        /* Else - We are at the first height, so use that for both above and
                  below */
    }
    /* Else - We are at the last height, so use that for both above and
              below */

    below_parameters[0] = modtran_results[below][LST_TRANSMISSION];
    below_parameters[1] = modtran_results[below][LST_UPWELLED_RADIANCE];
    below_parameters[2] = modtran_results[below][LST_DOWNWELLED_RADIANCE];

    if (above == below)
    {
        /* Use the below parameters since the same */
        at_height[0] = below_parameters[0];
        at_height[1] = below_parameters[1];
        at_height[2] = below_parameters[2];
    }
    else
    {
        /* Interpolate between the heights for each parameter */
        above_height = modtran_results[above][LST_HEIGHT];
        inv_height_diff = 1.0 / (above_height
                                 - modtran_results[below][LST_HEIGHT]);

        above_parameters[0] = modtran_results[above][LST_TRANSMISSION];
        above_parameters[1] = modtran_results[above][LST_UPWELLED_RADIANCE];
        above_parameters[2] = modtran_results[above][LST_DOWNWELLED_RADIANCE];

        for (i = 0; i < NUM_PARAMETERS; i++)
        {
            slope = (above_parameters[i] - below_parameters[i])
                    * inv_height_diff;

            intercept = above_parameters[i] - slope * above_height;

            at_height[i] = slope * interpolate_to + intercept;
        }
    }
}


/******************************************************************************
METHOD:  interpolate_to_location

PURPOSE: Interpolate to location of current pixel

******************************************************************************/
void interpolate_to_location
(
    REANALYSIS_POINTS *points,  /* I: The coordinate points */
    int *cell_vertices,         /* I: The vertices in the points to use */
    float **at_height,          /* I: current height atmospheric results */
    float interpolate_easting,  /* I: interpolate to easting */
    float interpolate_northing, /* I: interpolate to northing */
    float *parameters     /*O: interpolated pixel atmospheric parameters */
)
{
    int i, j;
    double inv_h[NUM_CELL_POINTS];
    double w[NUM_CELL_POINTS];
    double total = 0.0;

    /* shepard's method */
    for (i = 0; i < NUM_CELL_POINTS; i++)
    {
        inv_h[i] = 1.0 / sqrt (((points->utm_easting[cell_vertices[i]]
                                 - interpolate_easting)
                                * (points->utm_easting[cell_vertices[i]]
                                   - interpolate_easting))
                               +
                               ((points->utm_northing[cell_vertices[i]]
                                 - interpolate_northing)
                                * (points->utm_northing[cell_vertices[i]]
                                   - interpolate_northing)));

        total += inv_h[i];
    }

    /* Determine the weights for each vertex */
    for (i = 0; i < NUM_CELL_POINTS; i++)
    {
        w[i] = inv_h[i] / total;
    }

    /* For each parameter apply each vertex's weighted value */
    for (i = 0; i < NUM_PARAMETERS; i++)
    {
        parameters[i] = 0.0;
        for (j = 0; j < NUM_CELL_POINTS; j++)
        {
            parameters[i] += (w[j] * at_height[j][i]);
        }
    }
}


/*****************************************************************************
METHOD:  point_is_left_of_line

PURPOSE: Determines if a point is on the left side of the line or otherwise on
         the line or on  the right side of the line.

NOTE: This is based on a 2D geometry and when we are in UTM, that is the case.

RETURN: type = bool
    Value  Description
    -----  -------------------------------------------------------------------
    True   Indicates the value is on the left side of the line.
    False  Indicates the value is on the line or on the right side of the line.
*****************************************************************************/
bool point_is_left_of_line(int x0, int y0, int x1, int y1, int px, int py)
{
    float result = ((x1 - x0) * (py - y0)) - ((px - x0) * (y1 - y0));

    if (result > 0.0)
        return true;

    return false;
}


/*****************************************************************************
METHOD:  calculate_pixel_atmospheric_parameters

PURPOSE: Generate transmission, upwelled radiance, and downwelled radiance at
         each Landsat pixel

RETURN: SUCCESS
        FAILURE

*****************************************************************************/
int calculate_pixel_atmospheric_parameters
(
    Input_t *input,            /* I: input structure */
    REANALYSIS_POINTS *points, /* I: The coordinate points */
    char *xml_filename,        /* I: XML filename */
    char *dem_filename,        /* I: input DEM filename */
    char *emi_filename,        /* I: input Emissivity filename */
    float **modtran_results,   /* I: results from MODTRAN runs */
    bool verbose               /* I: value to indicate if intermediate
                                     messages be printed */
)
{
    char FUNC_NAME[] = "calculate_pixel_atmospheric_parameters";

    int line;
    int sample;
    int status;
    int closest_point;
    int offset;                 /* offset in the raw binary DEM file to seek to
                                   to begin reading the window in the DEM */
    bool first_sample;
    float current_easting;
    float current_northing;
    DISTANCE_ITEM *distances = NULL;

    int point;
    int vertex;
    int current_index;
    int cell_vertices[NUM_CELL_POINTS];

    float **at_height = NULL;
    float parameters[NUM_PARAMETERS];

    char *tmp_char = NULL;
    char scene_name[PATH_MAX];
    char thermal_filename[PATH_MAX];
    char upwelled_filename[PATH_MAX];
    char downwelled_filename[PATH_MAX];
    char transmittance_filename[PATH_MAX];

    FILE *dem_fd = NULL;
    FILE *thermal_fd = NULL;
    FILE *transmittance_fd = NULL;
    FILE *upwelled_fd = NULL;
    FILE *downwelled_fd = NULL;

    int16_t *dem = NULL;        /* input DEM data in meters */
    float *thermal_data;
    float **landsat_results;
    double current_height;
    char msg[MAX_STR_LEN];
    char *lst_data_dir = NULL;

    /* Use local variables for cleaner code */
    int num_cols = points->num_cols;
    int num_points = points->num_points;

    /* Grab the environment path to the LST_DATA_DIR */
    lst_data_dir = getenv ("LST_DATA_DIR");
    if (lst_data_dir == NULL)
    {
        RETURN_ERROR ("LST_DATA_DIR environment variable is not set",
                      FUNC_NAME, FAILURE);
    }

    /* Open the DEM for reading raw binary */
    dem_fd = fopen (dem_filename, "rb");
    if (dem_fd == NULL)
    {
        RETURN_ERROR ("Error opening the DEM file", FUNC_NAME, FAILURE);
    }

    /* Allocate memory for one line of the DEM */
    dem = calloc (input->thermal.size.s, sizeof (int16_t));
    if (dem == NULL)
    {
        RETURN_ERROR ("Error allocating memory for the DEM data",
                      FUNC_NAME, FAILURE);
    }

    /* Allocate memory for one line of the Thermal data */
    thermal_data = calloc (input->thermal.size.s, sizeof (float));
    if (thermal_data == NULL)
    {
        RETURN_ERROR ("Error allocating memory for the Thermal data",
                      FUNC_NAME, FAILURE);
    }

    /* Figure out the filenames */
    snprintf (scene_name, sizeof (scene_name), input->thermal.filename);
    tmp_char = strchr (scene_name, '_');
    if (tmp_char != NULL)
        *tmp_char = '\0';

    snprintf (thermal_filename, sizeof (thermal_filename),
              "%s_%s.img", scene_name, LST_THERMAL_RADIANCE_PRODUCT_NAME);
    snprintf (upwelled_filename, sizeof (upwelled_filename),
              "%s_%s.img", scene_name, LST_UPWELLED_RADIANCE_PRODUCT_NAME);
    snprintf (downwelled_filename, sizeof (downwelled_filename),
              "%s_%s.img", scene_name, LST_DOWNWELLED_RADIANCE_PRODUCT_NAME);
    snprintf (transmittance_filename, sizeof (transmittance_filename),
              "%s_%s.img", scene_name, LST_ATMOS_TRANS_PRODUCT_NAME);

    /* Open the intermediate binary files for writing
       Note: Needs to be deleted before release */
    thermal_fd = fopen (thermal_filename, "wb");
    if (thermal_fd == NULL)
    {
        sprintf (msg, "Opening report file: %s", thermal_filename);
        RETURN_ERROR (msg, FUNC_NAME, FAILURE);
    }

    transmittance_fd = fopen (transmittance_filename, "wb");
    if (transmittance_fd == NULL)
    {
        sprintf (msg, "Opening report file: %s", transmittance_filename);
        RETURN_ERROR (msg, FUNC_NAME, FAILURE);
    }

    upwelled_fd = fopen (upwelled_filename, "wb");
    if (upwelled_fd == NULL)
    {
        sprintf (msg, "Opening report file: %s", upwelled_filename);
        RETURN_ERROR (msg, FUNC_NAME, FAILURE);
    }

    downwelled_fd = fopen (downwelled_filename, "wb");
    if (downwelled_fd == NULL)
    {
        sprintf (msg, "Opening report file: %s", downwelled_filename);
        RETURN_ERROR (msg, FUNC_NAME, FAILURE);
    }

    /* Allocate memory for landsat_results */
    landsat_results = (float **) allocate_2d_array (3, input->thermal.size.s,
                                                    sizeof (float));
    if (landsat_results == NULL)
    {
        RETURN_ERROR ("Allocating landsat_results memory", FUNC_NAME, FAILURE);
    }

    /* Allocate memory for at_height */
    at_height = (float **) allocate_2d_array (NUM_CELL_POINTS, NUM_PARAMETERS,
                                              sizeof (float));
    if (at_height == NULL)
    {
        RETURN_ERROR ("Allocating at_height memory", FUNC_NAME, FAILURE);
    }

    /* Allocate memory to hold the distances to the first sample of data for
       the current line */
    distances = malloc (num_points * sizeof (DISTANCE_ITEM));
    if (distances == NULL)
    {
        ERROR_MESSAGE ("Allocating distances memory",
                       FUNC_NAME);
    }

    if (verbose)
    {
        LOG_MESSAGE ("Iterate through all lines in landsat scene\n",
                     FUNC_NAME);
    }

    /* Loop through each line in the image */
    for (line = 0; line < input->thermal.size.l; line++)
    {
        /* Print status on every 1000 lines */
        if (!(line % 1000))
        {
            if (verbose)
            {
                printf ("Processing line %d\r", line);
                fflush (stdout);
            }
        }

        /* Read the input thermal band data */
        if (!GetInputThermLine (input, line, thermal_data))
        {
            sprintf (msg, "Reading input thermal data for line %d", line);
            RETURN_ERROR (msg, FUNC_NAME, FAILURE);
        }

        /* Can also read in one line of DEM data here */
        /* Start reading DEM from the start_line */
        offset = sizeof (int16_t) * line * input->thermal.size.s;
        fseek (dem_fd, offset, SEEK_SET);
        if (fread (dem, sizeof (int16_t), input->thermal.size.s, dem_fd)
            != input->thermal.size.s)
        {
            sprintf (msg, "Error reading values from the DEM file "
                             "starting at line %d.", line);
            RETURN_ERROR (msg, FUNC_NAME, FAILURE);
        }

        /* Set first_sample to be true */
        first_sample = true;
        for (sample = 0; sample < input->thermal.size.s; sample++)
        {
            if (thermal_data[sample] != LST_NO_DATA_VALUE)
            {
                /* determine UTM coordinates of current pixel */
                current_easting = input->meta.ul_map_corner.x
                    + (sample * input->thermal.pixel_size[0]);
                current_northing = input->meta.ul_map_corner.y
                    - (line * input->thermal.pixel_size[1]);

                if (first_sample)
                {
                    /* compute distance between current pixel and each narr
                       point in UTM coordinates

                       Note: consider only calculating points within a small
                       range nearby */
                    for (point = 0; point < num_points; point++)
                    {
                        distance_in_utm (points->utm_easting[point],
                                         points->utm_northing[point],
                                         current_easting,
                                         current_northing,
                                         &distances[point].distance);

                        distances[point].index = point;
                    }

                    /* find the closest point */
                    qsort (distances, num_points, sizeof (DISTANCE_ITEM),
                           qsort_distance_compare_function);
                    closest_point = distances[0].index;
#if 0
snprintf (msg, sizeof (msg), "closest = %d", closest_point);
LOG_MESSAGE (msg, FUNC_NAME);
#endif

                    /* Now determine where we are in the point data cells */
                    if (point_is_left_of_line(
                        points->utm_easting[closest_point],
                        points->utm_northing[closest_point],
                        points->utm_easting[closest_point + num_cols],
                        points->utm_northing[closest_point + num_cols],
                        current_easting,
                        current_northing))
                    {
                        if (! point_is_left_of_line(
                            points->utm_easting[closest_point],
                            points->utm_northing[closest_point],
                            points->utm_easting[closest_point - 1],
                            points->utm_northing[closest_point - 1],
                            current_easting,
                            current_northing))
                        {
                            /* in quadrant top-left */
                            closest_point--;
                        }
                        else if (! point_is_left_of_line(
                            points->utm_easting[closest_point],
                            points->utm_northing[closest_point],
                            points->utm_easting[closest_point - num_cols],
                            points->utm_northing[closest_point - num_cols],
                            current_easting,
                            current_northing))
                        {
                            /* in quadrant bottom-left */
                            closest_point -= (num_cols + 1);
                        }
                        else
                        {
                            /* in quadrant bottom-right */
                            closest_point -= num_cols;
                        }
                    }
                    else if (! point_is_left_of_line(
                        points->utm_easting[closest_point],
                        points->utm_northing[closest_point],
                        points->utm_easting[closest_point + 1],
                        points->utm_northing[closest_point + 1],
                        current_easting,
                        current_northing))
                    {
                        if (point_is_left_of_line(
                            points->utm_easting[closest_point],
                            points->utm_northing[closest_point],
                            points->utm_easting[closest_point - num_cols],
                            points->utm_northing[closest_point - num_cols],
                            current_easting,
                            current_northing))
                        {
                            /* in quadrant bottom-right */
                            closest_point -= num_cols;
                        }
                        else
                        {
                            /* in quadrant bottom-left */
                            closest_point -= (num_cols + 1);
                        }
                    }
                    /* else in quadrant top-right */

                    /* determine index of four points in order to pull from
                       MODTRAN results */
                    /* LL */
                    cell_vertices[LL_POINT] = closest_point;
                    /* UL */
                    cell_vertices[UL_POINT] =
                        cell_vertices[LL_POINT] + num_cols;
                    /* UR */
                    cell_vertices[UR_POINT] = cell_vertices[UL_POINT] + 1;
                    /* LR */
                    cell_vertices[LR_POINT] = cell_vertices[LL_POINT] + 1;

                    /* Set first_sample to be false */
                    first_sample = false;
                }
                else
                {
                    /* Make sure we are:
                       left of LR->UR
                       and
                       below UR->UL

                       If not we need to advance to a different group of 4
                       points. */

                    if (! point_is_left_of_line(
                        points->utm_easting[cell_vertices[UR_POINT]],
                        points->utm_northing[cell_vertices[UR_POINT]],
                        points->utm_easting[cell_vertices[UL_POINT]],
                        points->utm_northing[cell_vertices[UL_POINT]],
                        current_easting,
                        current_northing))
                    {
                        /* We are above the line
                           Adjust the points and then test the right edge */

                        /* LL */
                        cell_vertices[LL_POINT] = cell_vertices[UL_POINT];
                        /* UL */
                        cell_vertices[UL_POINT] =
                            cell_vertices[LL_POINT] + num_cols;
                        /* UR */
                        cell_vertices[UR_POINT] = cell_vertices[UL_POINT] + 1;
                        /* LR */
                        cell_vertices[LR_POINT] = cell_vertices[LL_POINT] + 1;

                        /* Now test the right edge */
                        if (! point_is_left_of_line(
                            points->utm_easting[cell_vertices[LR_POINT]],
                            points->utm_northing[cell_vertices[LR_POINT]],
                            points->utm_easting[cell_vertices[UR_POINT]],
                            points->utm_northing[cell_vertices[UR_POINT]],
                            current_easting,
                            current_northing))
                        {
                            /* We are right of the line
                               Adjust the points */

                            /* LL */
                            cell_vertices[LL_POINT] = cell_vertices[LR_POINT];
                            /* UL */
                            cell_vertices[UL_POINT] =
                                cell_vertices[LL_POINT] + num_cols;
                            /* UR */
                            cell_vertices[UR_POINT] =
                                cell_vertices[UL_POINT] + 1;
                            /* LR */
                            cell_vertices[LR_POINT] =
                                cell_vertices[LL_POINT] + 1;
                        }
                        /* ELSE WE ARE OK FOR NOW */
                    }
                    else if (! point_is_left_of_line(
                        points->utm_easting[cell_vertices[LR_POINT]],
                        points->utm_northing[cell_vertices[LR_POINT]],
                        points->utm_easting[cell_vertices[UR_POINT]],
                        points->utm_northing[cell_vertices[UR_POINT]],
                        current_easting,
                        current_northing))
                    {
                        /* We are right of the line
                           Adjust the points and then test the top edge */

                        /* LL */
                        cell_vertices[LL_POINT] = cell_vertices[LR_POINT];
                        /* UL */
                        cell_vertices[UL_POINT] =
                            cell_vertices[LL_POINT] + num_cols;
                        /* UR */
                        cell_vertices[UR_POINT] =
                            cell_vertices[UL_POINT] + 1;
                        /* LR */
                        cell_vertices[LR_POINT] = cell_vertices[LL_POINT] + 1;

                        if (! point_is_left_of_line(
                            points->utm_easting[cell_vertices[UR_POINT]],
                            points->utm_northing[cell_vertices[UR_POINT]],
                            points->utm_easting[cell_vertices[UL_POINT]],
                            points->utm_northing[cell_vertices[UL_POINT]],
                            current_easting,
                            current_northing))
                        {
                            /* LL */
                            cell_vertices[LL_POINT] = cell_vertices[UL_POINT];
                            /* UL */
                            cell_vertices[UL_POINT] =
                                cell_vertices[LL_POINT] + num_cols;
                            /* UR */
                            cell_vertices[UR_POINT] =
                                cell_vertices[UL_POINT] + 1;
                            /* LR */
                            cell_vertices[LR_POINT] =
                                cell_vertices[LL_POINT] + 1;
                        }
                        /* ELSE WE ARE OK FOR NOW */
                    }
                    /* ELSE WE ARE OK FOR NOW */
                } /* END - Not first sample */

#if 0
snprintf (msg, sizeof (msg),
          "[%d, %d] cell_vertices = LL[%d], UL[%d], UR[%d], LR[%d]",
          line, sample, cell_vertices[LL_POINT], cell_vertices[UL_POINT],
                        cell_vertices[UR_POINT], cell_vertices[LR_POINT]);
LOG_MESSAGE (msg, FUNC_NAME);
#endif

                /* convert height from m to km */
                current_height = (double) dem[sample] / 1000.0;

                /* interpolate three parameters to that height at each of the
                   four closest points */

if (line == 3500 && sample == 4000)
{
    printf ("\n\n");
    printf ("easting [%lf]\n", current_easting);
    printf ("northing [%lf]\n", current_northing);
}
                for (vertex = 0; vertex < NUM_CELL_POINTS; vertex++)
                {
                    current_index = cell_vertices[vertex] * NUM_ELEVATIONS;

                    /* interpolate three atmospheric parameters to current
                       height */
                    interpolate_to_height (&modtran_results[current_index],
                                           current_height,
                                           at_height[vertex]);
//                                           &at_height[vertex][0]);
if (line == 3500 && sample == 4000)
{
    int bbbb;
    for (bbbb = 0; bbbb < NUM_ELEVATIONS; bbbb++)
    {
        printf ("%15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n",
                modtran_results[current_index+bbbb][LST_LATITUDE],
                modtran_results[current_index+bbbb][LST_LONGITUDE],
                modtran_results[current_index+bbbb][LST_HEIGHT],
                modtran_results[current_index+bbbb][LST_TRANSMISSION],
                modtran_results[current_index+bbbb][LST_UPWELLED_RADIANCE],
                modtran_results[current_index+bbbb][LST_DOWNWELLED_RADIANCE]);
    }

    printf ("at_height [%d][%d][%f] [%15.9f][%15.9f][%15.9f]\n",
            current_index, cell_vertices[vertex], current_height,
            at_height[vertex][0],
            at_height[vertex][1],
            at_height[vertex][2]);
}
                }

                /* interpolate parameters at appropriate height to location of
                   current pixel */
                interpolate_to_location (points, cell_vertices, at_height,
                                         current_easting, current_northing,
                                         &parameters[0]);

                /* convert radiances to W*m^(-2)*sr(-1) */
                landsat_results[0][sample] = parameters[0];
                landsat_results[1][sample] = parameters[1] * 10000.0;
                landsat_results[2][sample] = parameters[2] * 10000.0;

if (line == 3500 && sample == 4000)
{
    printf ("thermal [%f]\n", thermal_data[sample]);
    printf ("transmittance [%f]\n", landsat_results[0][sample]);
    printf ("transmittance [%f]\n", parameters[0]);
    printf ("upwelled [%f]\n", landsat_results[1][sample]);
    printf ("upwelled [%f]\n", parameters[1]);
    printf ("downwelled [%f]\n", landsat_results[2][sample]);
    printf ("downwelled [%f]\n", parameters[2]);
    printf ("\n\n");
}
            } /* END - if not FILL */
            else
            {
                landsat_results[0][sample] = LST_NO_DATA_VALUE;
                landsat_results[1][sample] = LST_NO_DATA_VALUE;
                landsat_results[2][sample] = LST_NO_DATA_VALUE;
            }
        } /* END - for sample */

        /* Write out the temporary binary output files
           Note: They need to be deleted before release */
        status = fwrite (thermal_data, sizeof (float),
                         input->thermal.size.s, thermal_fd);
        if (status != input->thermal.size.s)
        {
            sprintf (msg, "Writing to %s", thermal_filename);
            ERROR_MESSAGE (msg, FUNC_NAME);
        }

        status = fwrite (&landsat_results[0][0], sizeof (float),
                         input->thermal.size.s, transmittance_fd);
        if (status != input->thermal.size.s)
        {
            sprintf (msg, "Writing to %s", transmittance_filename);
            ERROR_MESSAGE (msg, FUNC_NAME);
        }

        status = fwrite (&landsat_results[1][0], sizeof (float),
                         input->thermal.size.s, upwelled_fd);
        if (status != input->thermal.size.s)
        {
            sprintf (msg, "Writing to %s", upwelled_filename);
            ERROR_MESSAGE (msg, FUNC_NAME);
        }

        status = fwrite (&landsat_results[2][0], sizeof (float),
                         input->thermal.size.s, downwelled_fd);
        if (status != input->thermal.size.s)
        {
            sprintf (msg, "Writing to %s", downwelled_filename);
            ERROR_MESSAGE (msg, FUNC_NAME);
        }
    } /* END - for line */

    /* Free allocated memory */
    free (distances);
    free (dem);
    free (thermal_data);
    thermal_data = NULL;

    status = free_2d_array ((void **) at_height);
    if (status != SUCCESS)
    {
        ERROR_MESSAGE ("Freeing memory: at_height\n", FUNC_NAME);
    }

    status = free_2d_array ((void **) landsat_results);
    if (status != SUCCESS)
    {
        ERROR_MESSAGE ("Freeing memory: landsat_results\n", FUNC_NAME);
    }

    /* Close the intermediate binary files
       Note: needs to be deleted before release */
    status = fclose (thermal_fd);
    if (status)
    {
        sprintf (msg, "Closing file %s", thermal_filename);
        ERROR_MESSAGE (msg, FUNC_NAME);
    }

    status = fclose (transmittance_fd);
    if (status)
    {
        sprintf (msg, "Closing file %s", transmittance_filename);
        ERROR_MESSAGE (msg, FUNC_NAME);
    }

    status = fclose (upwelled_fd);
    if (status)
    {
        sprintf (msg, "Closing file %s", upwelled_filename);
        ERROR_MESSAGE (msg, FUNC_NAME);
    }

    status = fclose (downwelled_fd);
    if (status)
    {
        sprintf (msg, "Closing file %s", downwelled_filename);
        ERROR_MESSAGE (msg, FUNC_NAME);
    }

    if (add_lst_band_product(xml_filename,
                             input->thermal.band_name,
                             LST_THERMAL_RADIANCE_PRODUCT_NAME,
                             LST_THERMAL_RADIANCE_BAND_NAME,
                             LST_THERMAL_RADIANCE_SHORT_NAME,
                             LST_THERMAL_RADIANCE_LONG_NAME,
                             LST_RADIANCE_UNITS,
                             0.0, 0.0) != SUCCESS)
    {
        ERROR_MESSAGE ("Failed adding LST band product", FUNC_NAME);
    }

    if (add_lst_band_product(xml_filename,
                             input->thermal.band_name,
                             LST_ATMOS_TRANS_PRODUCT_NAME,
                             LST_ATMOS_TRANS_BAND_NAME,
                             LST_ATMOS_TRANS_SHORT_NAME,
                             LST_ATMOS_TRANS_LONG_NAME,
                             LST_RADIANCE_UNITS,
                             0.0, 0.0) != SUCCESS)
    {
        ERROR_MESSAGE ("Failed adding LST band product", FUNC_NAME);
    }

    if (add_lst_band_product(xml_filename,
                             input->thermal.band_name,
                             LST_UPWELLED_RADIANCE_PRODUCT_NAME,
                             LST_UPWELLED_RADIANCE_BAND_NAME,
                             LST_UPWELLED_RADIANCE_SHORT_NAME,
                             LST_UPWELLED_RADIANCE_LONG_NAME,
                             LST_RADIANCE_UNITS,
                             0.0, 0.0) != SUCCESS)
    {
        ERROR_MESSAGE ("Failed adding LST band product", FUNC_NAME);
    }

    if (add_lst_band_product(xml_filename,
                             input->thermal.band_name,
                             LST_DOWNWELLED_RADIANCE_PRODUCT_NAME,
                             LST_DOWNWELLED_RADIANCE_BAND_NAME,
                             LST_DOWNWELLED_RADIANCE_SHORT_NAME,
                             LST_DOWNWELLED_RADIANCE_LONG_NAME,
                             LST_RADIANCE_UNITS,
                             0.0, 0.0) != SUCCESS)
    {
        ERROR_MESSAGE ("Failed adding LST band product", FUNC_NAME);
    }

    return SUCCESS;
}
