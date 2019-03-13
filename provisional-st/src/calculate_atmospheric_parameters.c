#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <getopt.h>
#include <errno.h>
#include <string.h>
#include <float.h>


#include "const.h"
#include "utilities.h"
#include "input.h"
#include "output.h"
#include "intermediate_data.h"
#include "grid_points.h"
#include "integrate.h"
#include "interpolate.h"
#include "modtran_utils.h"
#include "calculate_atmospheric_parameters.h"

/*****************************************************************************
DESCRIPTION: Using values produced by MODTRAN runs at grid points, calculate
atmospheric transmission, upwelled radiance, and downwelled radiance for each
pixel in the Landsat scene.  Also, create bands with these values.  A thermal
radiance band is also created based on a Landsat thermal band and parameters.
*****************************************************************************/


/*****************************************************************************
METHOD:  planck_eq

PURPOSE: Using Planck's equation to calculate radiance at each wavelength for
         current temperature.
*****************************************************************************/
static void planck_eq
(
    double *wavelength, /* I: Each wavelength */
    int num_elements,   /* I: Number of wavelengths to calculate */
    double temperature, /* I: The temperature to calculate for */
    double *bb_radiance /* O: the blackbody results for each wavelength */
)
{
    int i;

    /* Planck Const hecht pg, 585 ## units: Js */
    double PLANCK_CONST = 6.6260755e-34;

    /* Boltzmann Gas Const halliday et 2001 -- units: J/K */
    double BOLTZMANN_GAS_CONST = 1.3806503e-23;

    /* Speed of Light -- units: um/s */
    double SPEED_OF_LIGHT = 299792458e6;
    double SPEED_OF_LIGHT_SQRD = SPEED_OF_LIGHT * SPEED_OF_LIGHT;

    for (i = 0; i < num_elements; i++)
    {
        /* Compute the Planck Blackbody Eq [W/m^2 sr um].
           Convert to W/cm^2 sr micron to match modtran units
           (additional 1e-4 multiplier). */
        bb_radiance[i] = 2e8*PLANCK_CONST*SPEED_OF_LIGHT_SQRD
                       * pow(wavelength[i], -5.0)
                       / (exp(PLANCK_CONST*SPEED_OF_LIGHT
                              /(wavelength[i]*BOLTZMANN_GAS_CONST*temperature))
                          - 1.0);
    }
}


/* Set up structure for managing the blackbody radiation work arrays.*/
struct rad_work_array
{
    int size;
    double *rad;
    double *prod;
};

static struct rad_work_array rad_work = {0, NULL, NULL};
    
static void free_rad_workspace()
{
    free(rad_work.rad);
    free(rad_work.prod);
    rad_work.rad = NULL;
    rad_work.prod = NULL;
    rad_work.size = 0;
}


/*****************************************************************************
MODULE:  calculate_lt

PURPOSE: Calculate blackbody radiance from temperature using spectral response
         function.

RETURN: SUCCESS
        FAILURE
*****************************************************************************/
static int calculate_lt
(
    double temperature,         /*I: temperature */
    double **spectral_response, /*I: spectral response function */
    int num_srs,                /*I: number of spectral response points */
    double *radiance            /*O: blackbody radiance */
)
{
    char FUNC_NAME[] = "calculate_lt";
    int i;
    double rs_integral;
    double temp_integral;
    double *product;
    double *blackbody_radiance;

     /* Allocate memory */
    if (num_srs > rad_work.size)
    {
        rad_work.rad = realloc(rad_work.rad, num_srs*sizeof(double));
        rad_work.prod = realloc(rad_work.prod, num_srs*sizeof(double));
        if (rad_work.rad == NULL || rad_work.prod == NULL)
        {
            free_rad_workspace();
            RETURN_ERROR ("Allocating blackbody workspace memory", FUNC_NAME,
                          FAILURE);
        }
        rad_work.size = num_srs;
    }
    blackbody_radiance = rad_work.rad;
    product = rad_work.prod;

    /* integrate spectral response over wavelength */
    if (int_tabulated (spectral_response[0], spectral_response[1], num_srs,
                       &rs_integral) != SUCCESS)
    {
        RETURN_ERROR ("Calling int_tabulated\n", FUNC_NAME, FAILURE);
    }

    /* Use planck's blackbody radiance equation to calculate radiance at each
       wavelength for the current temperature */
    planck_eq (spectral_response[0], num_srs, temperature, blackbody_radiance);

    /* Multiply the calculated planck radiance by the spectral response and
       integrate over wavelength to get one number for current temp */
    for (i = 0; i < num_srs; i++)
    {
        product[i] = blackbody_radiance[i] * spectral_response[1][i];
    }

    if (int_tabulated (spectral_response[0], product, num_srs,
                       &temp_integral) != SUCCESS)
    {
        RETURN_ERROR ("Calling int_tabulated\n", FUNC_NAME, FAILURE);
    }

    /* Divide above result by integral of spectral response function */
    *radiance = temp_integral / rs_integral;

    return SUCCESS;
}


/*****************************************************************************
MODULE:  linear_interpolate_over_modtran

PURPOSE: Simulate IDL (interpol) function for ST.
*****************************************************************************/
static void linear_interpolate_over_modtran
(
    double (*modtran)[4], /* I: The MODTRAN data - provides both the a and b */
    int index,        /* I: The MODTRAN temperatur to use for a */
    double *c,        /* I: The Landsat wavelength grid points */
    int num_in,       /* I: Number of input data and grid points*/
    int num_out,      /* I: Number of output grid points */
    double *x         /* O: Interpolated output results */
)
{
    int i;
    int o;

    double d1 = 0.0;
    double d2 = 0.0;
    double g;
    double g1 = 0.0;
    double g2 = 0.0;

    int a = index; /* MODTRAN radiance for specififc temp */
    int b = 0;     /* MODTRAN wavelength */

    for (o = 0; o < num_out; o++)
    {
        g = c[o];

        /* Initialize to the first two */
        d1 = modtran[0][a];
        d2 = modtran[1][a];
        g1 = modtran[0][b];
        g2 = modtran[1][b];

        for (i = 0; i < num_in-1; i++)
        {
            if (g <= modtran[i][b] && g > modtran[i+1][b])
            {
                /* Found it in the middle of the data */
                d1 = modtran[i][a];
                d2 = modtran[i+1][a];
                g1 = modtran[i][b];
                g2 = modtran[i+1][b];
                break;
            }
        }

        if (i == num_in-1)
        {
            /* Less than the last so use the last two */
            d1 = modtran[i-1][a];
            d2 = modtran[i][a];
            g1 = modtran[i-1][b];
            g2 = modtran[i][b];
        }

        /* Apply the formula for linear interpolation */
        x[o] = d1 + (g - g1) / (g2 - g1) * (d2 - d1);
    }
}


/*****************************************************************************
MODULE:  calculate_lobs

PURPOSE: Calculate observed radiance from MODTRAN results and the spectral
         response function.

RETURN: SUCCESS
        FAILURE
*****************************************************************************/
static int calculate_lobs
(
    double (*modtran)[4],       /*I: MODTRAN results with wavelengths */
    double **spectral_response, /*I: spectral response function */
    int num_entries,            /*I: number of MODTRAN points */
    int num_srs,                /*I: number of spectral response points */
    int index,                  /*I: column index for data be used */
    double *radiance            /*O: LOB outputs */
)
{
    char FUNC_NAME[] = "calculate_lobs";
    int i;
    double *temp_rad;
    double rs_integral;
    double temp_integral;
    double *product;

    /* Allocate memory */
    if (num_srs > rad_work.size)
    {
        rad_work.rad = realloc(rad_work.rad, num_srs*sizeof(double));
        rad_work.prod = realloc(rad_work.prod, num_srs*sizeof(double));
        if (rad_work.rad == NULL || rad_work.prod == NULL)
        {
            free_rad_workspace();
            RETURN_ERROR ("Allocating blackbody workspace memory", FUNC_NAME,
                          FAILURE);
        }
        rad_work.size = num_srs;
    }
    temp_rad = rad_work.rad;
    product = rad_work.prod;

    /* Integrate spectral response over wavelength */
    if (int_tabulated (spectral_response[0], spectral_response[1], num_srs,
                       &rs_integral) != SUCCESS)
    {
        RETURN_ERROR ("Calling int_tabulated\n", FUNC_NAME, FAILURE);
    }

    /* Interpolate MODTRAN radiance to Landsat wavelengths */
    linear_interpolate_over_modtran (modtran, index, spectral_response[0],
                                     num_entries, num_srs, temp_rad);

    /* Multiply the calculated radiance by the spectral response and integrate
       over wavelength to get one number for current temperature */
    for (i = 0; i < num_srs; i++)
    {
        product[i] = temp_rad[i] * spectral_response[1][i];
    }

    if (int_tabulated (spectral_response[0], product, num_srs,
                       &temp_integral) != SUCCESS)
    {
        RETURN_ERROR ("Calling int_tabulated\n", FUNC_NAME, FAILURE);
    }

    /* Divide above result by integral of spectral response function */
    *radiance = temp_integral / rs_integral;

    return SUCCESS;
}


/*****************************************************************************
METHOD:  calculate_point_atmospheric_parameters

PURPOSE: Generate transmission, upwelled radiance, and downwelled radiance at
         each height for each NARR point that is used.

RETURN: SUCCESS
        FAILURE
*****************************************************************************/
static int calculate_point_atmospheric_parameters
(
    Input_Data_t *input,       /* I: Input structure */
    GRID_POINTS *grid_points,  /* I: The coordinate points */
    MODTRAN_POINTS *modtran_results /* I/O: Atmospheric parameters from 
                                   MODTRAN */
)
{
    char FUNC_NAME[] = "calculate_point_atmospheric_parameters";

    FILE *fd;
    FILE *used_points_fd;

    int i;
    int j;
    int entry;

    double *spectral_response[2];
    double temp_radiance_0;
    double obs_radiance_0;
    double temp_radiance_273;
    double temp_radiance_310;
    double delta_radiance_inv; /* inverse of radiance differences; used to
                                  compute transmittance and upwelled radiance */
    int counter;
    int index;
    int num_entries;   /* Number of MODTRAN output results to read and use */
    int num_srs;       /* Number of spectral response values available */

    char *st_data_dir = NULL;
    char current_file[PATH_MAX]; /* Used for MODTRAN info (input), MODTRAN data
                          (input), and atmospheric parameters (output) files */
    char srs_file_path[PATH_MAX];
    char msg[PATH_MAX];

    double modtran_wavelength;
    double modtran_radiance;
    double zero_temp;
    double (*current_data)[4] = NULL;
    int max_radiance_record_count = 0; /* max number of radiance records read */
    double y_0;
    double y_1;
    double tau; /* Transmission */
    double lu;  /* Upwelled Radiance */
    double ld;  /* Downwelled Radiance */

    /* Temperature and albedo */
    int temperature[3] = { 273, 310, 000 };
    double albedo[3] = { 0.0, 0.0, 0.1 };


    st_data_dir = getenv ("ST_DATA_DIR");
    if (st_data_dir == NULL)
    {
        RETURN_ERROR ("ST_DATA_DIR environment variable is not set",
                      FUNC_NAME, FAILURE);
    }

    /* Allocate memory for maximum spectral response count */
    spectral_response[0] = malloc(MAX_SRS_COUNT*sizeof(double));
    spectral_response[1] = malloc(MAX_SRS_COUNT*sizeof(double));
    if (spectral_response[0] == NULL || spectral_response[1] == NULL)
    {
        RETURN_ERROR ("Allocating spectral_response memory",
                      FUNC_NAME, FAILURE);
    }

    /* Determine the spectral response file to read */
    if (input->meta.instrument == INST_TM
        && input->meta.satellite == SAT_LANDSAT_4)
    {
        num_srs = L4_TM_SRS_COUNT;

        snprintf (srs_file_path, sizeof (srs_file_path),
                  "%s/%s", st_data_dir, "L4_Spectral_Response.txt");
    }
    else if (input->meta.instrument == INST_TM
        && input->meta.satellite == SAT_LANDSAT_5)
    {
        num_srs = L5_TM_SRS_COUNT;

        snprintf (srs_file_path, sizeof (srs_file_path),
                  "%s/%s", st_data_dir, "L5_Spectral_Response.txt");
    }
    else if (input->meta.instrument == INST_ETM
             && input->meta.satellite == SAT_LANDSAT_7)
    {
        num_srs = L7_TM_SRS_COUNT;

        snprintf (srs_file_path, sizeof (srs_file_path),
                  "%s/%s", st_data_dir, "L7_Spectral_Response.txt");
    }
    else if (input->meta.instrument == INST_OLI_TIRS
             && input->meta.satellite == SAT_LANDSAT_8)
    {
        num_srs = L8_OLITIRS_SRS_COUNT;

        snprintf (srs_file_path, sizeof (srs_file_path),
                  "%s/%s", st_data_dir, "L8_Spectral_Response.txt");
    }
    else
    {
        RETURN_ERROR ("invalid instrument type", FUNC_NAME, FAILURE);
    }

    /* Read the selected spectral response file */
    snprintf (msg, sizeof (msg),
              "Reading Spectral Response File [%s]", srs_file_path);
    LOG_MESSAGE (msg, FUNC_NAME);
    fd = fopen (srs_file_path, "r");
    if (fd == NULL)
    {
        RETURN_ERROR ("Can't open Spectral Response file", FUNC_NAME, FAILURE);
    }

    for (i = 0; i < num_srs; i++)
    {
        if (fscanf (fd, "%lf %lf%*c", &spectral_response[0][i],
                    &spectral_response[1][i]) == EOF)
        {
            RETURN_ERROR ("Failed reading spectral response file",
                          FUNC_NAME, FAILURE);
        }
    }
    fclose (fd);

    /* Calculate Lt for each specific temperature */
    if (calculate_lt (273, spectral_response, num_srs, &temp_radiance_273)
        != SUCCESS)
    {
        RETURN_ERROR ("Calling calculate_lt for 273K", FUNC_NAME, FAILURE);
    }
    if (calculate_lt (310, spectral_response, num_srs, &temp_radiance_310)
        != SUCCESS)
    {
        RETURN_ERROR ("Calling calculate_lt for 310K", FUNC_NAME, FAILURE);
    }

    /* Compute the multiplier for the transmittance and upwelled radiance
       calculations in the following loop. */
    delta_radiance_inv = 1/(temp_radiance_310 - temp_radiance_273);

    /* Output information about the used points, primarily useful for
       plotting them against the scene */
    used_points_fd = fopen ("used_points.txt", "w");
    if (used_points_fd == NULL)
    {
        RETURN_ERROR ("Can't open used_points.txt file",
                      FUNC_NAME, FAILURE);
    }

    /* Iterate through all grid points and heights */
    counter = 0;
    for (i = 0; i < grid_points->count; i++)
    {
        GRID_POINT *grid_point = &grid_points->points[i];
        MODTRAN_POINT *modtran_point = & modtran_results->points[i];

        /* Don't process the points that didn't have a MODTRAN run. */
        if (!modtran_point->ran_modtran)
        {
            continue;
        }

        fprintf (used_points_fd, "\"%d\"|\"%f\"|\"%f\"\n",
                 i, grid_point->map_x, grid_point->map_y);

        for (j = 0; j < modtran_point->count; j++)
        {
            /* Read the st_modtran.info file for the 000 execution
               (when MODTRAN is run at 0K)
               We read the zero_temp from this file, and also the record count
               The record count is the same for all three associated runs */
            snprintf(current_file, sizeof(current_file), 
                "%03d_%03d_%03d_%03d/%1.3f/000/0.1/st_modtran.hdr",
                grid_point->row, grid_point->col, grid_point->narr_row,
                grid_point->narr_col,
                modtran_point->elevations[j].elevation_directory);

            fd = fopen (current_file, "r");
            if (fd == NULL)
            {
                snprintf (msg, sizeof (msg),
                          "Can't open MODTRAN information file [%s]",
                           current_file);
                RETURN_ERROR (msg, FUNC_NAME, FAILURE);
            }
            /* Retrieve the temperature from this lowest atmospheric layer */
            if (fscanf (fd, "%*s %lf%*c", &zero_temp) != 1)
            {
                RETURN_ERROR ("End of file (EOF) is met before"
                              " reading TARGET_PIXEL_SURFACE_TEMPERATURE",
                              FUNC_NAME, FAILURE);
            }
            /* Determine number of entries in current file */
            if (fscanf (fd, "%*s %d%*c", &num_entries) != 1)
            {
                RETURN_ERROR ("End of file (EOF) is met before"
                              " reading RADIANCE_RECORD_COUNT",
                              FUNC_NAME, FAILURE);
            }
            fclose (fd);

            /* For each height, read in radiance information for three
               MODTRAN runs.  Columns of array are organized as follows:
               wavelength | 273,0.0 | 310,0.0 | 000,0.1 */
            if (num_entries > max_radiance_record_count)
            {
                max_radiance_record_count = num_entries;
                current_data = realloc(current_data,
                                       num_entries*sizeof(double[4]));
                if (current_data == NULL)
                {
                    RETURN_ERROR ("Allocating current_data memory",
                                  FUNC_NAME, FAILURE);
                }
            }

            /* Iterate through the three pairs of parameters */
            for (index = 1; index < 4; index++)
            {
                /* Define MODTRAN data file */
                snprintf(current_file, sizeof(current_file), 
                    "%03d_%03d_%03d_%03d/%1.3f/%03d/%1.1f/st_modtran.data",
                    grid_point->row, grid_point->col, grid_point->narr_row,
                    grid_point->narr_col,
                    modtran_point->elevations[j].elevation_directory,
                    temperature[index - 1],
                    albedo[index - 1]);

                fd = fopen (current_file, "r");
                if (fd == NULL)
                {
                    RETURN_ERROR ("Can't open MODTRAN data file",
                                  FUNC_NAME, FAILURE);
                }
                for (entry = 0; entry < num_entries; entry++)
                {
                    if (fscanf (fd, "%lf %lf%*c",
                                &modtran_wavelength, &modtran_radiance)
                        != 2)
                    {
                        RETURN_ERROR ("Failed reading st_modtran.dat lines",
                                      FUNC_NAME, FAILURE);
                    }

                    /* If we are on the first file set the wavelength value
                       for the data array */
                    if (index == 1)
                    {
                        current_data[entry][0] = modtran_wavelength;
                    }
                    /* Place radiance into data array for current point at
                       current height */
                    current_data[entry][index] = modtran_radiance;

                }
                fclose (fd);

                counter++;
            }

            /* Parameters from 3 MODTRAN runs
               Lobs = Lt*tau + Lu; m = tau; b = Lu; */
            if (calculate_lobs (current_data, spectral_response,
                                num_entries, num_srs, 1, &y_0)
                != SUCCESS)
            {
                RETURN_ERROR ("Calling calculate_lobs for height y_0",
                              FUNC_NAME, FAILURE);
            }

            if (calculate_lobs (current_data, spectral_response,
                                num_entries, num_srs, 2, &y_1)
                != SUCCESS)
            {
                RETURN_ERROR ("Calling calculate_lobs for height y_1",
                              FUNC_NAME, FAILURE);
            }

            tau = (y_1 - y_0)*delta_radiance_inv; /* Transmittance */
            lu = (temp_radiance_310*y_0 - temp_radiance_273*y_1)
               * delta_radiance_inv;  /* Upwelled Radiance */

            /* Determine Lobs and Lt when MODTRAN was run at 0K - calculate 
               downwelled */
            if (calculate_lt (zero_temp, spectral_response, num_srs,
                              &temp_radiance_0) != SUCCESS)
            {
                RETURN_ERROR ("Calling calculate_lt for zero temp (0Kelvin)",
                              FUNC_NAME, FAILURE);
            }

            if (calculate_lobs (current_data, spectral_response,
                                num_entries, num_srs, 3, &obs_radiance_0)
                != SUCCESS)
            {
                RETURN_ERROR ("Calling calculate_lobs for (0Kelvin)",
                              FUNC_NAME, FAILURE);
            }

            /* Calculate the downwelled radiance. These are all equivalent:
               Ld = (((Lobs - Lu) / tau)
                     - (Lt * WATER_EMISSIVITY)) / (1.0 - WATER_EMISSIVITY)
               Ld = (((Lobs - Lu) / tau)
                     - (Lt * WATER_EMISSIVITY)) / WATER_ALBEDO
               Ld = (((Lobs - Lu) / tau)
                     - (Lt * WATER_EMISSIVITY)) * INV_WATER_ALBEDO */
            ld = (((obs_radiance_0 - lu) / tau)
                  - (temp_radiance_0 * WATER_EMISSIVITY)) * INV_WATER_ALBEDO;

            /* Place results into MODTRAN results array */
            modtran_point->elevations[j].transmission = tau;
            modtran_point->elevations[j].upwelled_radiance = lu;
            modtran_point->elevations[j].downwelled_radiance = ld;
        } /* END - modtran_point->count loop */
    } /* END - count loop */
    fclose (used_points_fd);

    /* Free allocated memory */
    free(current_data);
    current_data = NULL;
    free(spectral_response[0]);
    spectral_response[0] = NULL;
    free(spectral_response[1]);
    spectral_response[1] = NULL;
    free_integration_workspace();
    free_rad_workspace();

    /* Write atmospheric transmission, upwelled radiance, and downwelled 
       radiance for each elevation for each point to a file */
    snprintf (current_file, sizeof (current_file),
              "atmospheric_parameters.txt");
    snprintf (msg, sizeof (msg),
              "Creating Atmospheric Parameters File = [%s]\n", current_file);
    LOG_MESSAGE (msg, FUNC_NAME);
    fd = fopen (current_file, "w");
    if (fd == NULL)
    {
        RETURN_ERROR ("Can't open atmospheric_parameters.txt file",
                      FUNC_NAME, FAILURE);
    }
    for (i = 0; i < grid_points->count; i++)
    {
        MODTRAN_POINT *modtran_point = & modtran_results->points[i];

        /* Only write parameters for grid points where MODTRAN was run */
        if (!modtran_point->ran_modtran)
        {
            continue;
        }

        for (j = 0; j < modtran_point->count; j++)
        {
            fprintf (fd, "%f,%f,%12.9f,%12.9f,%12.9f,%12.9f\n",
                 modtran_point->lat,
                 modtran_point->lon,
                 modtran_point->elevations[j].elevation,
                 modtran_point->elevations[j].transmission,
                 modtran_point->elevations[j].upwelled_radiance,
                 modtran_point->elevations[j].downwelled_radiance);
        }
    }
    fclose (fd);

    return SUCCESS;
}


/*****************************************************************************
METHOD:  calculate_pixel_atmospheric_parameters

PURPOSE: Generate transmission, upwelled radiance, and downwelled radiance at
         each Landsat pixel

RETURN: SUCCESS
        FAILURE
*****************************************************************************/
static int calculate_pixel_atmospheric_parameters
(
    Input_Data_t *input,       /* I: input structure */
    GRID_POINTS *points,       /* I: coordinate points */
    char *xml_filename,        /* I: XML filename */
    Espa_internal_meta_t xml_metadata, /* I: XML metadata */
    MODTRAN_POINTS *modtran_results /* I: results from MODTRAN runs */
)
{
    char FUNC_NAME[] = "calculate_pixel_atmospheric_parameters";

    int line;
    int sample;
    int ipt;                   /* grid point index */

    bool first_sample;

    double easting;
    double northing;

    Geoloc_t *space = NULL;    /* Geolocation information */
    Space_def_t space_def;     /* Space definition (projection values) */
    Img_coord_float_t img;     /* Floating point image coordinates */
    Geo_coord_t geo;           /* Geodetic coordinates */
    double *grid_pt_lon;       /* grid point longitudes (radians) */
    double *grid_pt_lat;       /* grid point latitudes (radians) */


    GRID_ITEM *grid_points = NULL;

    int center_point;
    int cell_vertices[NUM_CELL_POINTS];

    double parameters[AHP_NUM_PARAMETERS];
    double avg_distance_ll;
    double avg_distance_ul;
    double avg_distance_ur;
    double avg_distance_lr;

    Intermediate_Data_t inter;

    int16_t *elevation_data = NULL; /* input elevation data in meters */

    double current_height;
    char msg[MAX_STR_LEN];

    /* Use local variables for cleaner code */
    int num_cols = points->cols;
    int num_points = points->count;
    int pixel_count = input->lines * input->samples;
    int pixel_loc;

    /* Open the intermedate data files */
    if (open_intermediate(input, &inter) != SUCCESS)
    {
        RETURN_ERROR("Opening intermediate data files", FUNC_NAME, FAILURE);
    }

    /* Allocate memory for the intermedate data */
    if (allocate_intermediate(&inter, pixel_count) != SUCCESS)
    {
        RETURN_ERROR("Allocating memory for intermediate data",
                     FUNC_NAME, FAILURE);
    }

    /* Allocate memory for elevation */
    elevation_data = malloc(pixel_count*sizeof(int16_t));
    if (elevation_data == NULL)
    {
        RETURN_ERROR("Allocating elevation_data memory", FUNC_NAME, FAILURE);
    }

    /* Allocate memory to hold the grid_points to the first sample of data for
       the current line */
    grid_points = malloc (num_points * sizeof (GRID_ITEM));
    if (grid_points == NULL)
    {
        RETURN_ERROR ("Allocating grid_points memory", FUNC_NAME, FAILURE);
    }

    /* Convert the grid point latitudes and longitudes to radians for use
       in calcuations. */
    grid_pt_lon = malloc(num_points * sizeof(double));
    grid_pt_lat = malloc(num_points * sizeof(double));
    if (grid_pt_lon == NULL || grid_pt_lat == NULL)
        RETURN_ERROR ("Allocating grid point lat/long memory", FUNC_NAME,
                      FAILURE);
    for (ipt = 0; ipt < num_points; ipt++)
    {
        grid_pt_lon[ipt] = points->points[ipt].lon*RADIANS_PER_DEGREE;
        grid_pt_lat[ipt] = points->points[ipt].lat*RADIANS_PER_DEGREE;
    }

    /* Read thermal and elevation data into memory */
    if (read_input(input, inter.band_thermal, elevation_data, pixel_count)
        != SUCCESS)
    {
        RETURN_ERROR ("Reading thermal and elevation bands", FUNC_NAME,
                      FAILURE);
    }

    /* Get geolocation space definition */
    if (!get_geoloc_info(&xml_metadata, &space_def))
    {
        RETURN_ERROR ("Getting space metadata from XML file", FUNC_NAME,
                     FAILURE);
    }
    space = setup_mapping(&space_def);
    if (space == NULL)
    {
        RETURN_ERROR ("Setting up geolocation mapping", FUNC_NAME, FAILURE);
    }

    /* Show some status messages */
    LOG_MESSAGE("Iterate through all pixels in Landsat scene", FUNC_NAME);
    snprintf(msg, sizeof(msg), "Pixel Count = %d", pixel_count);
    LOG_MESSAGE(msg, FUNC_NAME);
    snprintf(msg,  sizeof(msg),"Lines = %d, Samples = %d",
             input->lines, input->samples);
    LOG_MESSAGE(msg, FUNC_NAME);

    /* Loop through each line in the image */
    for (line = 0, pixel_loc = 0; line < input->lines; line++)
    {
        /* Print status on every 1000 lines */
        if (!(line % 1000))
        {
            printf ("Processing line %d\n", line);
            fflush (stdout);
        }

        northing = input->meta.ul_map_corner.y - line*input->y_pixel_size;

        /* Set first_sample to be true */
        first_sample = true;
        for (sample = 0; sample < input->samples; sample++, pixel_loc++)
        {
            /* Skip fill pixels. */
            if (inter.band_thermal[pixel_loc] == ST_NO_DATA_VALUE)
            {
                inter.band_upwelled[pixel_loc] = ST_NO_DATA_VALUE;
                inter.band_downwelled[pixel_loc] = ST_NO_DATA_VALUE;
                inter.band_transmittance[pixel_loc] = ST_NO_DATA_VALUE;

#if OUTPUT_CELL_DESIGNATION_BAND
                inter.band_cell[pixel_loc] = 0;
#endif
                continue;
            }

            /* Determine latitude and longitude for current line/sample */
            img.l = line;
            img.s = sample;
            img.is_fill = false;
            if (!from_space(space, &img, &geo))
            {
                RETURN_ERROR ("Mapping from line/sample to longitude/"
                              "latitude", FUNC_NAME, FAILURE);
            }
            easting = input->meta.ul_map_corner.x + sample*input->x_pixel_size;

            if (first_sample)
            {
                /* Determine the first center point from all of the
                   available points */
                center_point = determine_first_center_grid_point(
                    grid_pt_lon, grid_pt_lat, geo.lon,
                    geo.lat, num_points, grid_points);

                /* Set first_sample to be false */
                first_sample = false;
            }
            else
            {
                /* Determine the center point from the current 9 grid
                   points for the current line/sample */
                center_point = determine_center_grid_point(
                    grid_pt_lon, grid_pt_lat, geo.lon,
                    geo.lat, NUM_GRID_POINTS, grid_points);
            }

            /* Fix the index values, since the points are from a new line
               or were messed up during determining the center point */
            grid_points[CC_GRID_POINT].index = center_point;
            grid_points[LL_GRID_POINT].index = center_point - 1 - num_cols;
            grid_points[LC_GRID_POINT].index = center_point - 1;
            grid_points[UL_GRID_POINT].index = center_point - 1 + num_cols;
            grid_points[UC_GRID_POINT].index = center_point + num_cols;
            grid_points[UR_GRID_POINT].index = center_point + 1 + num_cols;
            grid_points[RC_GRID_POINT].index = center_point + 1;
            grid_points[LR_GRID_POINT].index = center_point + 1 - num_cols;
            grid_points[DC_GRID_POINT].index = center_point - num_cols;

            /* Fix the distances, since the points are from a new line or
               were messed up during determining the center point */
            determine_grid_point_distances(grid_pt_lon, grid_pt_lat,
                                           geo.lon, geo.lat,
                                           NUM_GRID_POINTS, grid_points);

            /* Determine the average distances for each quadrant around
               the center point. We only need to use the three outer grid 
               points.
               The true average values are the values computed
               below divided by 3.  But since we're only using these these
               values comparatively, we can forego the division. */
            avg_distance_ll = grid_points[DC_GRID_POINT].distance
                            + grid_points[LL_GRID_POINT].distance
                            + grid_points[LC_GRID_POINT].distance;

            avg_distance_ul = grid_points[LC_GRID_POINT].distance
                            + grid_points[UL_GRID_POINT].distance
                            + grid_points[UC_GRID_POINT].distance;

            avg_distance_ur = grid_points[UC_GRID_POINT].distance
                            + grid_points[UR_GRID_POINT].distance
                            + grid_points[RC_GRID_POINT].distance;

            avg_distance_lr = grid_points[RC_GRID_POINT].distance
                            + grid_points[LR_GRID_POINT].distance
                            + grid_points[DC_GRID_POINT].distance;

            /* Determine which quadrant is closer and setup the cell
               vertices to interpolate over based on that */
            if (avg_distance_ll < avg_distance_ul
                && avg_distance_ll < avg_distance_ur
                && avg_distance_ll < avg_distance_lr)
            { /* LL Cell */
                cell_vertices[LL_POINT] = center_point - 1 - num_cols;
            }
            else if (avg_distance_ul < avg_distance_ur
                     && avg_distance_ul < avg_distance_lr)
            { /* UL Cell */
                cell_vertices[LL_POINT] = center_point - 1;
            }
            else if (avg_distance_ur < avg_distance_lr)
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
            inter.band_cell[pixel_loc] = cell_vertices[LL_POINT];
#endif

            /* Convert height from m to km -- Same as 1.0 / 1000.0 */
            current_height = (double) elevation_data[pixel_loc] * 0.001;

            /* Interpolate three parameters to the height and location of
               the current pixel. */
            interpolate_parameters(modtran_results->points, points,
                                   cell_vertices, current_height, easting,
                                   northing, &parameters[0]);

            /* Convert radiances to W*m^(-2)*sr(-1) */
            inter.band_upwelled[pixel_loc] =
                parameters[AHP_UPWELLED_RADIANCE] * 10000.0;
            inter.band_downwelled[pixel_loc] =
                parameters[AHP_DOWNWELLED_RADIANCE] * 10000.0;
            inter.band_transmittance[pixel_loc] =
                parameters[AHP_TRANSMISSION];
        } /* END - for sample */
    } /* END - for line */

    /* Write out the temporary intermediate output files */
    if (write_intermediate(&inter, pixel_count) != SUCCESS)
    {
        sprintf (msg, "Writing to intermediate data files");
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }

    /* Free allocated memory */
    free(grid_points);
    free(grid_pt_lat);
    free(grid_pt_lon);
    free(elevation_data);
    free_intermediate(&inter);

    /* Close the intermediate binary files */
    if (close_intermediate(&inter) != SUCCESS)
    {
        sprintf (msg, "Closing file intermediate data files");
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }

    /* Add the ST intermediate bands to the metadata file */
    if (add_st_band_product(xml_filename,
                             input->reference_band_name,
                             inter.thermal_filename,
                             ST_THERMAL_RADIANCE_PRODUCT_NAME,
                             ST_THERMAL_RADIANCE_BAND_NAME,
                             ST_THERMAL_RADIANCE_SHORT_NAME,
                             ST_THERMAL_RADIANCE_LONG_NAME,
                             ST_RADIANCE_UNITS,
                             0.0, 0.0) != SUCCESS)
    {
        ERROR_MESSAGE ("Failed adding ST thermal radiance band product", 
            FUNC_NAME);
    }

    if (add_st_band_product(xml_filename,
                             input->reference_band_name,
                             inter.transmittance_filename,
                             ST_ATMOS_TRANS_PRODUCT_NAME,
                             ST_ATMOS_TRANS_BAND_NAME,
                             ST_ATMOS_TRANS_SHORT_NAME,
                             ST_ATMOS_TRANS_LONG_NAME,
                             ST_RADIANCE_UNITS,
                             0.0, 0.0) != SUCCESS)
    {
        ERROR_MESSAGE ("Failed adding ST atmospheric transmission band "
            "product", FUNC_NAME);
    }

    if (add_st_band_product(xml_filename,
                             input->reference_band_name,
                             inter.upwelled_filename,
                             ST_UPWELLED_RADIANCE_PRODUCT_NAME,
                             ST_UPWELLED_RADIANCE_BAND_NAME,
                             ST_UPWELLED_RADIANCE_SHORT_NAME,
                             ST_UPWELLED_RADIANCE_LONG_NAME,
                             ST_RADIANCE_UNITS,
                             0.0, 0.0) != SUCCESS)
    {
        ERROR_MESSAGE ("Failed adding ST upwelled radiance band product", 
            FUNC_NAME);
    }

    if (add_st_band_product(xml_filename,
                             input->reference_band_name,
                             inter.downwelled_filename,
                             ST_DOWNWELLED_RADIANCE_PRODUCT_NAME,
                             ST_DOWNWELLED_RADIANCE_BAND_NAME,
                             ST_DOWNWELLED_RADIANCE_SHORT_NAME,
                             ST_DOWNWELLED_RADIANCE_LONG_NAME,
                             ST_RADIANCE_UNITS,
                             0.0, 0.0) != SUCCESS)
    {
        ERROR_MESSAGE ("Failed adding ST downwelled radiance band product", 
            FUNC_NAME);
    }

    return SUCCESS;
}


/****************************************************************************
Method: usage

Description: Display help/usage information to the user.
****************************************************************************/
void usage()
{
    printf("Surface Temperature - st_atmospheric_parameters\n");
    printf("\n");
    printf("Generates interpolated atmospheric parameters covering the scene"
           " data.\n");
    printf("\n");
    printf("usage: st_atmospheric_parameters"
           " --xml=<filename>"
           " [--debug]\n");
    printf("\n");
    printf ("where the following parameters are required:\n");
    printf ("    --xml: name of the input XML file\n");
    printf ("\n");
    printf ("where the following parameters are optional:\n");
    printf ("    --debug: should debug output be generated?"
            " (default is false)\n");
    printf ("\n");
    printf ("st_atmospheric_parameters --help will print the ");
    printf ("usage statement\n");
    printf ("\n");
    printf ("Example: st_atmospheric_parameters"
            " --xml=LE07_L1T_028031_20041227_20160513_01_T1.xml\n");
    printf ("Note: This application must run from the directory"
            " where the input data is located.\n\n");
}


/*****************************************************************************
Method:  get_args

Description:  Gets the command-line arguments and validates that the required
              arguments were specified.

Returns: Type = int
    Value           Description
    -----           -----------
    FAILURE         Error getting the command-line arguments or a command-line
                    argument and associated value were not specified
    SUCCESS         No errors encountered

Notes:
    1. Memory is allocated for the input and output files.  All of these
       should be character pointers set to NULL on input.  The caller is
       responsible for freeing the allocated memory upon successful return.
*****************************************************************************/
static int get_args
(
    int argc,           /* I: number of cmd-line args */
    char *argv[],       /* I: string of cmd-line args */
    char *xml_filename, /* I: address of input XML metadata filename  */
    bool *debug         /* O: debug flag */
)
{
    int c;                         /* current argument index */
    int option_index;              /* index of the command line option */
    static int debug_flag = 0;     /* debug flag */
    char errmsg[MAX_STR_LEN];      /* error message */
    char FUNC_NAME[] = "get_args"; /* function name */

    static struct option long_options[] = {
        {"debug", no_argument, &debug_flag, 1},
        {"xml", required_argument, 0, 'i'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    /* Loop through all the cmd-line options */
    opterr = 0; /* turn off getopt_long error msgs as we'll print our own */
    while (1)
    {
        /* optstring in call to getopt_long is empty since we will only
           support the long options */
        c = getopt_long(argc, argv, "", long_options, &option_index);
        if (c == -1)
        {
            /* Out of cmd-line options */
            break;
        }

        switch (c)
        {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;

            case 'h':              /* help */
                usage ();
                return FAILURE;
                break;

            case 'i':              /* xml infile */
                snprintf(xml_filename, PATH_MAX, "%s", optarg);
                break;

            case '?':
            default:
                snprintf(errmsg, sizeof(errmsg),
                         "Unknown option %s", argv[optind - 1]);
                usage();
                RETURN_ERROR(errmsg, FUNC_NAME, FAILURE);
                break;
        }
    }

    /* Make sure the XML file was specified */
    if (strlen(xml_filename) <= 0)
    {
        usage();
        RETURN_ERROR("XML input file is a required argument", FUNC_NAME,
                     FAILURE);
    }

    /* Set the debug flag */
    if (debug_flag)
        *debug = true;
    else
        *debug = false;

    return SUCCESS;
}


/*****************************************************************************
Method:  main

Description:  Main for the application.
*****************************************************************************/
int main(int argc, char *argv[])
{
    char FUNC_NAME[] = "main";

    Espa_internal_meta_t xml_metadata;  /* XML metadata structure */
    char xml_filename[PATH_MAX];        /* Input XML filename */
    bool debug;                         /* Debug flag for debug output */
    Input_Data_t *input = NULL;         /* Input data and meta data */
    GRID_POINTS grid_points;            /* NARR grid points */
    MODTRAN_POINTS modtran_points;      /* Points that are processed through
                                           MODTRAN */

    /* Read the command-line arguments */
    if (get_args(argc, argv, xml_filename, &debug)
        != SUCCESS)
    {
        RETURN_ERROR("calling get_args", FUNC_NAME, EXIT_FAILURE);
    }

    /* Validate the input metadata file */
    if (validate_xml_file(xml_filename) != SUCCESS)
    {
        /* Error messages already written */
        return EXIT_FAILURE;
    }

    /* Initialize the metadata structure */
    init_metadata_struct(&xml_metadata);

    /* Parse the metadata file into our internal metadata structure; also
       allocates space as needed for various pointers in the global and band
       metadata */
    if (parse_metadata(xml_filename, &xml_metadata) != SUCCESS)
    {
        /* Error messages already written */
        return EXIT_FAILURE;
    }

    /* Open input file, read metadata, and set up buffers */
    input = open_input(&xml_metadata);
    if (input == NULL)
    {
        RETURN_ERROR("opening input files", FUNC_NAME, EXIT_FAILURE);
    }

    /* Load the grid points */
    if (load_grid_points(&grid_points) != SUCCESS)
    {
        RETURN_ERROR("calling load_grid_points", FUNC_NAME, EXIT_FAILURE);
    }

    /* Allocate and initialize the memory need to hold the MODTRAN results */
    if (initialize_modtran_points(&grid_points, &modtran_points) != SUCCESS)
    {
        RETURN_ERROR("calling initializing_modtran_points", FUNC_NAME, 
            EXIT_FAILURE);
    }

    /* Generate parameters for each height and NARR point */
    if (calculate_point_atmospheric_parameters(input, &grid_points, 
        &modtran_points) != SUCCESS)
    {
        RETURN_ERROR("calling calculate_point_atmospheric_parameters",
            FUNC_NAME, EXIT_FAILURE);
    }

    /* Process the grid points */
    printf("%d %d %d\n", grid_points.count, grid_points.rows, grid_points.cols);
    printf("%d %d %d\n", grid_points.points[0].index, grid_points.points[0].row,
        grid_points.points[0].col);
    printf("%d %d %d\n", grid_points.points[1].index, grid_points.points[1].row,
        grid_points.points[1].col);
    printf("%d %d %d\n", grid_points.points[2].index, grid_points.points[2].row,
        grid_points.points[2].col);
    printf("%d %d %d\n", grid_points.points[3].index, grid_points.points[3].row,
        grid_points.points[3].col);
    printf("%d %d %d\n", grid_points.points[15].index, 
       grid_points.points[15].row, grid_points.points[15].col);

    /* Using the values made at the grid points, generate atmospheric 
       parameters for each Landsat pixel */ 
    if (calculate_pixel_atmospheric_parameters(input, &grid_points, 
        xml_filename, xml_metadata, &modtran_points) != SUCCESS)
    {
        RETURN_ERROR("calling calculate_pixel_atmospheric_parameters", 
            FUNC_NAME, EXIT_FAILURE);
    }

    /* Free metadata */
    free_metadata(&xml_metadata);

    /* Free the grid and MODTRAN points */
    free_grid_points(&grid_points);
    free_modtran_points(&modtran_points);

    /* Close the input file and free the structure */
    close_input(input);

    return EXIT_SUCCESS;
}
