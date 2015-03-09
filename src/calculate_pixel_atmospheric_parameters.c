
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
#define NUM_INTERMEDIATE_BANDS 3


#define OUTPUT_CELL_DESIGNATION_BAND 1
#define OUTPUT_INTERMEDIATE_BANDS 1


/* Defines the distance to the current pixel, along with the index of the
   point
   So that we can find the index of the closest point to start determining the
   correct cell to use */
typedef struct
{
    int index;
    double distance;
} GRID_ITEM;


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
} GRID_POINTS;


/* A qsort routine that can be used with the GRID_ITEM items */
int qsort_grid_compare_function
(
    const void *grid_item_a,
    const void *grid_item_b
)
{
    double a = (*(GRID_ITEM*)grid_item_a).distance;
    double b = (*(GRID_ITEM*)grid_item_b).distance;

    if (a < b)
        return -1;
    else if (b < a)
        return 1;

    return 0;
}


/******************************************************************************
METHOD:  distance_in_utm

PURPOSE: Calculate distances between UTM coordiantes

RETURN: double - The distance.

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
double distance_in_utm
(
    double e0,
    double n0,
    double e2,
    double n2
)
{
    /* The UTM coordinates we are using have the 500000 false easting applied
       to them, so we need to remove that before applying the distance
       calculation. */
    double e0_adj;
    double e1_term;
    double e2_adj;

    double sr_e0;
    double sr_e1;
    double sr_e2;

    double edist;

    e0_adj = e0 - UTM_FALSE_EASTING;
    e2_adj = e2 - UTM_FALSE_EASTING;
    e1_term = (e0_adj + e2_adj) * INV_TWO;

    sr_e0 = UTM_SCALE_FACTOR / (cos (e0_adj * INV_UTM_EQUATORIAL_RADIUS));

    sr_e1 = UTM_SCALE_FACTOR / (cos (e1_term * INV_UTM_EQUATORIAL_RADIUS));

    sr_e2 = UTM_SCALE_FACTOR / (cos (e2_adj * INV_UTM_EQUATORIAL_RADIUS));

    edist = ((e2 - e0) * INV_SIX)
            * (sr_e0 + 4.0 * sr_e1 + sr_e2);

    return sqrt (edist * edist + (n2 - n0) * (n2 - n0));
}


/******************************************************************************
METHOD:  interpolate_to_height

PURPOSE: Interpolate to height of current pixel

******************************************************************************/
void interpolate_to_height
(
    double **modtran_results, /* I: results from MODTRAN runs for a point */
    double interpolate_to,    /* I: current landsat pixel height */
    double *at_height         /* O: interpolated height for point */
)
{
    int parameter;
    int elevation;
    int below = 0;
    int above = 0;

    double below_parameters[NUM_PARAMETERS];
    double above_parameters[NUM_PARAMETERS];

    double slope;
    double intercept;

    double above_height;
    double inv_height_diff; /* To remove the multiple divisions */

    /* Find the height to use that is below the interpolate_to height */
    for (elevation = 0; elevation < NUM_ELEVATIONS; elevation++)
    {
        if (modtran_results[elevation][LST_HEIGHT] < interpolate_to)
        {
            below = elevation; /* Last match will always be the one we want */
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

        for (parameter = 0; parameter < NUM_PARAMETERS; parameter++)
        {
            slope = (above_parameters[parameter] - below_parameters[parameter])
                    * inv_height_diff;

            intercept = above_parameters[parameter] - slope * above_height;

            at_height[parameter] = slope * interpolate_to + intercept;
        }
    }
}


/******************************************************************************
METHOD:  interpolate_to_location

PURPOSE: Interpolate to location of current pixel

******************************************************************************/
void interpolate_to_location
(
    REANALYSIS_POINTS *points,   /* I: The coordinate points */
    int *cell_vertices,          /* I: The vertices in the points to use */
    double **at_height,          /* I: current height atmospheric results */
    double interpolate_easting,  /* I: interpolate to easting */
    double interpolate_northing, /* I: interpolate to northing */
    double *parameters     /*O: interpolated pixel atmospheric parameters */
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
    double result = ((x1 - x0) * (py - y0)) - ((px - x0) * (y1 - y0));

    if (result > 0.0)
        return true;

    return false;
}


/*****************************************************************************
METHOD:  determine_grid_point_distances

PURPOSE: Determines the distances for the current set of grid points.

NOTE: The indexes of the grid points are assumed to be populated.

RETURN: None
*****************************************************************************/
void determine_grid_point_distances
(
    REANALYSIS_POINTS *points, /* I: All the available points */
    double easting,            /* I: Easting of the current line/sample */
    double northing,           /* I: Northing of the current line/sample */
    int num_grid_points,       /* I: The number of grid points to operate on */
    GRID_ITEM *grid_points     /* I/O: Sorted to determine the center grid
                                       point */
)
{
    int point;

    /* Populate the distances to the grid points */
    for (point = 0; point < num_grid_points; point++)
    {
        grid_points[point].distance = distance_in_utm (
            points->utm_easting[grid_points[point].index],
            points->utm_northing[grid_points[point].index],
            easting, northing);
    }
}


/*****************************************************************************
METHOD:  determine_center_grid_point

PURPOSE: Determines the index of the center point from the current set of grid
         points.

NOTE: The indexes of the grid points are assumed to be populated.

RETURN: type = int
    Value  Description
    -----  -------------------------------------------------------------------
    index  The index of the center point.
*****************************************************************************/
int determine_center_grid_point
(
    REANALYSIS_POINTS *points, /* I: All the available points */
    double easting,            /* I: Easting of the current line/sample */
    double northing,           /* I: Northing of the current line/sample */
    int num_grid_points,       /* I: The number of grid points to operate on */
    GRID_ITEM *grid_points     /* I/O: Sorted to determine the center grid
                                       point */
)
{
    determine_grid_point_distances (points, easting, northing,
                                    num_grid_points, grid_points);

    /* Sort them to find the closest one */
    qsort (grid_points, num_grid_points, sizeof (GRID_ITEM),
           qsort_grid_compare_function);

    return grid_points[0].index;
}


/*****************************************************************************
METHOD:  determine_first_center_grid_point

PURPOSE: Determines the index of the first center point to use for the current
         line.  Only called when the fist valid point for a line is
         encountered.  The point is determined from all of the available
         points.

RETURN: type = int
    Value  Description
    -----  -------------------------------------------------------------------
    index  The index of the center point.
*****************************************************************************/
int determine_first_center_grid_point
(
    REANALYSIS_POINTS *points, /* I: All the available points */
    double easting,            /* I: Easting of the current line/sample */
    double northing,           /* I: Northing of the current line/sample */
    GRID_ITEM *grid_points     /* I/O: Memory passed in, polulated and
                                       sorted to determine the center grid
                                       point */
)
{
    int point;

    /* Assign the point indexes for all grid points */
    for (point = 0; point < points->num_points; point++)
    {
        grid_points[point].index = point;
    }

    return determine_center_grid_point (points, easting, northing,
                                        points->num_points, grid_points);
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
    double **modtran_results,  /* I: results from MODTRAN runs */
    bool verbose               /* I: value to indicate if intermediate
                                     messages be printed */
)
{
    char FUNC_NAME[] = "calculate_pixel_atmospheric_parameters";

    int line;
    int sample;
    int status;
    int offset;                 /* offset in the raw binary DEM file to seek to
                                   to begin reading the window in the DEM */
    bool first_sample;
    double easting;
    double northing;
    GRID_ITEM *grid_points = NULL;

    int vertex;
    int current_index;
    int center_point;
    int cell_vertices[NUM_CELL_POINTS];

    double **at_height = NULL;
    double parameters[NUM_PARAMETERS];
    double avg_distance_ll;
    double avg_distance_ul;
    double avg_distance_ur;
    double avg_distance_lr;

    char *tmp_char = NULL;
    char scene_name[PATH_MAX];

#if OUTPUT_INTERMEDIATE_BANDS
    char thermal_filename[PATH_MAX];
    char upwelled_filename[PATH_MAX];
    char downwelled_filename[PATH_MAX];
    char transmittance_filename[PATH_MAX];
    FILE *thermal_fd = NULL;
    FILE *transmittance_fd = NULL;
    FILE *upwelled_fd = NULL;
    FILE *downwelled_fd = NULL;
#endif

    FILE *dem_fd = NULL;
    int16_t *dem = NULL;        /* input DEM data in meters */
    float *thermal_data;
    float **landsat_results;

#if OUTPUT_CELL_DESIGNATION_BAND
    char cell_filename[PATH_MAX];
    FILE *cell_fd = NULL;
    uint8_t *cell_data;
#endif

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

#if OUTPUT_CELL_DESIGNATION_BAND
    /* Allocate memory for one line of the CELL data */
    cell_data = calloc (input->thermal.size.s, sizeof (uint8_t));
    if (cell_data == NULL)
    {
        RETURN_ERROR ("Error allocating memory for the CELL data",
                      FUNC_NAME, FAILURE);
    }
#endif

    /* Figure out the scene name to apply to the filenames */
    snprintf (scene_name, sizeof (scene_name), input->thermal.filename);
    tmp_char = strchr (scene_name, '_');
    if (tmp_char != NULL)
        *tmp_char = '\0';

#if OUTPUT_CELL_DESIGNATION_BAND
    snprintf (cell_filename, sizeof (cell_filename),
              "%s_%s.img", scene_name, "cells");
#endif

#if OUTPUT_INTERMEDIATE_BANDS
    snprintf (thermal_filename, sizeof (thermal_filename),
              "%s_%s.img", scene_name, LST_THERMAL_RADIANCE_PRODUCT_NAME);
    snprintf (upwelled_filename, sizeof (upwelled_filename),
              "%s_%s.img", scene_name, LST_UPWELLED_RADIANCE_PRODUCT_NAME);
    snprintf (downwelled_filename, sizeof (downwelled_filename),
              "%s_%s.img", scene_name, LST_DOWNWELLED_RADIANCE_PRODUCT_NAME);
    snprintf (transmittance_filename, sizeof (transmittance_filename),
              "%s_%s.img", scene_name, LST_ATMOS_TRANS_PRODUCT_NAME);
#endif

    /* Open the intermediate binary files for writing */
#if OUTPUT_CELL_DESIGNATION_BAND
    cell_fd = fopen (cell_filename, "wb");
    if (cell_fd == NULL)
    {
        sprintf (msg, "Opening report file: %s", thermal_filename);
        RETURN_ERROR (msg, FUNC_NAME, FAILURE);
    }
#endif

#if OUTPUT_INTERMEDIATE_BANDS
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
#endif

    /* Allocate memory for landsat_results */
    landsat_results = (float **) allocate_2d_array (3, input->thermal.size.s,
                                                    sizeof (float));
    if (landsat_results == NULL)
    {
        RETURN_ERROR ("Allocating landsat_results memory", FUNC_NAME, FAILURE);
    }

    /* Allocate memory for at_height */
    at_height = (double **) allocate_2d_array (NUM_CELL_POINTS, NUM_PARAMETERS,
                                               sizeof (double));
    if (at_height == NULL)
    {
        RETURN_ERROR ("Allocating at_height memory", FUNC_NAME, FAILURE);
    }

    /* Allocate memory to hold the grid_points to the first sample of data for
       the current line */
    grid_points = malloc (num_points * sizeof (GRID_ITEM));
    if (grid_points == NULL)
    {
        ERROR_MESSAGE ("Allocating grid_points memory", FUNC_NAME);
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
                /* Determine UTM coordinates for current line/sample */
                easting = input->meta.ul_map_corner.x
                    + (sample * input->thermal.pixel_size[0]);
                northing = input->meta.ul_map_corner.y
                    - (line * input->thermal.pixel_size[1]);

                if (first_sample)
                {
                    /* Determine the first center point from all of the
                       available points */
                    center_point = determine_first_center_grid_point(
                                       points, easting, northing, grid_points);

                    /* Set first_sample to be false */
                    first_sample = false;
                }
                else
                {
                    /* Determine the center point from the current 9 grid
                       points for the current line/sample */
                    center_point = determine_center_grid_point(
                                       points, easting, northing,
                                       NUM_GRID_POINTS, grid_points);
                }

                /* Fix the index values, since the points are from a new line
                   or were messed up during determining the center point */
                grid_points[CC_GRID_POINT].index = center_point;
                grid_points[LL_GRID_POINT].index =
                    center_point - 1 - num_cols;
                grid_points[LC_GRID_POINT].index = center_point - 1;
                grid_points[UL_GRID_POINT].index =
                    center_point - 1 + num_cols;
                grid_points[UC_GRID_POINT].index =
                    center_point + num_cols;
                grid_points[UR_GRID_POINT].index =
                    center_point + 1 + num_cols;
                grid_points[RC_GRID_POINT].index = center_point + 1;
                grid_points[LR_GRID_POINT].index =
                    center_point + 1 - num_cols;
                grid_points[DC_GRID_POINT].index =
                    center_point - num_cols;

                /* Fix the distances, since the points are from a new line or
                   were messed up during determining the center point */
                determine_grid_point_distances (points, easting, northing,
                                                NUM_GRID_POINTS, grid_points);

                /* Determine the average distances for each quadrant around
                   the center point
                   We only need to use the three outer grid points */
                avg_distance_ll = (grid_points[DC_GRID_POINT].distance
                                   + grid_points[LL_GRID_POINT].distance
                                   + grid_points[LC_GRID_POINT].distance)
                                  / 3.0;

                avg_distance_ul = (grid_points[LC_GRID_POINT].distance
                                   + grid_points[UL_GRID_POINT].distance
                                   + grid_points[UC_GRID_POINT].distance)
                                  / 3.0;

                avg_distance_ur = (grid_points[UC_GRID_POINT].distance
                                   + grid_points[UR_GRID_POINT].distance
                                   + grid_points[RC_GRID_POINT].distance)
                                  / 3.0;

                avg_distance_lr = (grid_points[RC_GRID_POINT].distance
                                   + grid_points[LR_GRID_POINT].distance
                                   + grid_points[DC_GRID_POINT].distance)
                                  / 3.0;

                /* Determine which quadrant is closer and setup the cell
                   vertices to interpolate over based on that */
                if (avg_distance_ll < avg_distance_ul
                    && avg_distance_ll < avg_distance_ur
                    && avg_distance_ll < avg_distance_lr)
                { /* LL Cell */
                    cell_vertices[LL_POINT] = center_point - 1 - num_cols;
                }
                else if (avg_distance_ul < avg_distance_ll
                    && avg_distance_ul < avg_distance_ur
                    && avg_distance_ul < avg_distance_lr)
                { /* UL Cell */
                    cell_vertices[LL_POINT] = center_point - 1;
                }
                else if (avg_distance_ur < avg_distance_ll
                    && avg_distance_ur < avg_distance_ul
                    && avg_distance_ur < avg_distance_lr)
                { /* UR Cell */
                    cell_vertices[LL_POINT] = center_point;
                }
                else
                { /* LR Cell */
                    cell_vertices[LL_POINT] = center_point - num_cols;
                }

                /* UL Point */
                cell_vertices[UL_POINT] = cell_vertices[LL_POINT] + num_cols;
                /* UR Point */
                cell_vertices[UR_POINT] = cell_vertices[UL_POINT] + 1;
                /* LR Point */
                cell_vertices[LR_POINT] = cell_vertices[LL_POINT] + 1;

#if OUTPUT_CELL_DESIGNATION_BAND
                cell_data[sample] = cell_vertices[LL_POINT];
#endif

                /* convert height from m to km -- Same as 1.0 / 1000.0 */
                current_height = (double) dem[sample] * 0.001;

                /* interpolate three parameters to that height at each of the
                   four closest points */
                for (vertex = 0; vertex < NUM_CELL_POINTS; vertex++)
                {
                    current_index = cell_vertices[vertex] * NUM_ELEVATIONS;

                    /* interpolate three atmospheric parameters to current
                       height */
                    interpolate_to_height (&modtran_results[current_index],
                                           current_height,
                                           at_height[vertex]);
                }

                /* interpolate parameters at appropriate height to location of
                   current pixel */
                interpolate_to_location (points, cell_vertices, at_height,
                                         easting, northing, &parameters[0]);

                /* convert radiances to W*m^(-2)*sr(-1) */
                landsat_results[0][sample] = parameters[0];
                landsat_results[1][sample] = parameters[1] * 10000.0;
                landsat_results[2][sample] = parameters[2] * 10000.0;
            } /* END - if not FILL */
            else
            {
#if OUTPUT_CELL_DESIGNATION_BAND
                cell_data[sample] = 0;
#endif
                landsat_results[0][sample] = LST_NO_DATA_VALUE;
                landsat_results[1][sample] = LST_NO_DATA_VALUE;
                landsat_results[2][sample] = LST_NO_DATA_VALUE;
            }
        } /* END - for sample */

        /* Write out the temporary binary output files */
#if OUTPUT_CELL_DESIGNATION_BAND
        status = fwrite (cell_data, sizeof (uint8_t),
                         input->thermal.size.s, cell_fd);
        if (status != input->thermal.size.s)
        {
            sprintf (msg, "Writing to %s", thermal_filename);
            ERROR_MESSAGE (msg, FUNC_NAME);
        }
#endif

#if OUTPUT_INTERMEDIATE_BANDS
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
#endif
    } /* END - for line */

    /* Free allocated memory */
    free (grid_points);
    free (dem);
#if OUTPUT_CELL_DESIGNATION_BAND
    free (cell_data);
#endif
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

    /* Close the intermediate binary files */
#if OUTPUT_CELL_DESIGNATION_BAND
    status = fclose (cell_fd);
    if (status)
    {
        sprintf (msg, "Closing file %s", cell_filename);
        ERROR_MESSAGE (msg, FUNC_NAME);
    }
#endif

#if OUTPUT_INTERMEDIATE_BANDS
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
#endif

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
