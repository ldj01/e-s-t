#include <stdarg.h>
#include <math.h>


#include "const.h"
#include "2d_array.h"
#include "utilities.h"
#include "input.h"
#include "lst_types.h"


/******************************************************************************
MODULE:  planck_eq

PURPOSE: using planck's equaiton to calculate radiance at each wavelength for
         current temperature

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/30/2014   Song Guo         Original Development
******************************************************************************/
int planck_eq
(
    float *wavelength,
    int num_srs,
    float temperature,
    float *black_radiance
)
{
    char FUNC_NAME[] = "planck_eq";
    int i;
    float *lambda;
    float h;
    float k;
    float c;

    /* Allocate memory */
    lambda = (float *) malloc (num_srs * sizeof (float));
    if (lambda == NULL)
    {
        RETURN_ERROR ("Allocating lambda memory", FUNC_NAME, FAILURE);
    }

    /* Planck Const hecht pg, 585 ## units: Js */
    h = 6.6260755 * (float) 10e-34;

    /* Boltzmann Gas Const halliday et 2001 ## units: J/K */
    k = 1.3806503 * (float) 10e-23;

    /* Speed of Light ## units: m/s */
    c = 299792458.0;

    for (i = 0; i < num_srs; i++)
    {
        /* Lambda intervals of Landsat5 spectral response locations microns 
           units: m */
        lambda[i] = wavelength[i] * (float) 10e-6;

        /* Compute the Planck Blackbody Eq [W/m^2 sr um] */
        black_radiance[i] = 2.0 * h * pow (c, 2.0)
                            * ((float) (10e-6) * pow (lambda[i], -5.0))
                            * (1.0 / (exp ((h * c)
                                            / (lambda [i] * k * temperature))
                                      - 1.0));

        /* convert to W/cm^2 sr micron to match modtran units */
        black_radiance[i] /= 10000.0;
    }
    free (lambda);

    return SUCCESS;
}


/******************************************************************************
MODULE:  spline

PURPOSE: spline constructs a cubic spline given a set of x and y values,
         through these values.

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/29/2014   Song Guo         Modified from Numerical Recipes in C
                             (ISBN 0-521-43108-5)
******************************************************************************/
int spline
(
    float *x,
    float *y,
    int n,
    float yp1,
    float ypn,
    float *y2
)
{
    char FUNC_NAME[] = "spline";
    int i;
    double p;
    double qn;
    double sig;
    double un;
    double *u = NULL;

    u = (double *) malloc ((unsigned) (n - 1) * sizeof (double));
    if (u == NULL)
    {
        RETURN_ERROR ("Can't allocate memory", FUNC_NAME, FAILURE);
    }

    /* Set the lower boundary */
    if (yp1 > 0.99e30)
    {
        /* To be "natural" */
        y2[0] = 0.0;
        u[0] = 0.0;
    }
    else
    {
        /* To have a specified first derivative */
        y2[0] = -0.5;
        u[0] = (3.0 / (x[1] - x[0]))
               * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
    }

    /* Set the upper boundary */
    if (ypn > 0.99e30)
    {
        /* To be "natural" */
        qn = 0.0;
        un = 0.0;
    }
    else
    {
        /* To have a specified first derivative */
        qn = 0.5;
        un = (3.0 / (x[n - 1] - x[n - 2]))
             * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
    }

    /* Perform decomposition of the tridiagonal algorithm */
    for (i = 1; i <= n - 2; i++)
    {
        sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);

        p = sig * y2[i - 1] + 2.0;

        y2[i] = (sig - 1.0) / p;

        u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
               - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);

        u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
    }
    y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);

    /* Perform the backsubstitution of the tridiagonal algorithm */
    for (i = n - 2; i >= 0; i--)
    {
        y2[i] = y2[i] * y2[i + 1] + u[i];
    }

    free (u);

    return SUCCESS;
}


/******************************************************************************
MODULE:  splint

PURPOSE: splint uses the cubic spline generated with spline to interpolate
         values in the XY table

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/29/2014   Song Guo         Modified from online code
******************************************************************************/
static int klo = -1;
static int khi = -1;
static double one_sixth = (1.0 / 6.0); /* To remove a division */
void splint
(
    float *xa,
    float *ya,
    float *y2a,
    int n,
    float x,
    float *y
)
{
    int k;
    double h;
    double b;
    double a;

    if (klo < 0)
    {
        klo = 0;
        khi = n - 1;
    }
    else
    {
        if (x < xa[klo])
            klo = 0;
        if (x > xa[khi])
            khi = n - 1;
    }

    while (khi - klo > 1)
    {
        k = (khi + klo) >> 1;

        if (xa[k] > x)
            khi = k;
        else
            klo = k;
    }

    h = xa[khi] - xa[klo];

    if (h == 0.0)
    {
        *y = 0.0;
    }
    else
    {
        a = (xa[khi] - x) / h;

        b = (x - xa[klo]) / h;

        *y = a * ya[klo]
             + b * ya[khi]
             + ((a * a * a - a) * y2a[klo]
                + (b * b * b - b) * y2a[khi]) * (h * h) * one_sixth;
    }
}


/******************************************************************************
MODULE:  int_tabulated

PURPOSE: This function integrates a tabulated set of data { x(i) , f(i) },
         on the closed interval [min(X) , max(X)].

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/29/2014   Song Guo         Original Development
******************************************************************************/
int int_tabulated
(
    float *x,         /*I: The tabulated X-value data */
    float *f,         /*I: The tabulated F-value data */
    int nums,         /*I: number of points */
    float *result_out /*O: integraeted result */
)
{
    char FUNC_NAME[] = "int_tabulated";
    float *temp = NULL;
    float *z = NULL;
    float xmin;
    float xmax;
    int i;
    int *ii = NULL;
    int ii_count;
    float h;
    float result;
    int segments;

    /* Allocate memory */
    temp = (float *) malloc (nums * sizeof (float));
    if (temp == NULL)
    {
        RETURN_ERROR ("Allocating temp memory", FUNC_NAME, FAILURE);
    }

    z = (float *) malloc (nums * sizeof (float));
    if (z == NULL)
    {
        RETURN_ERROR ("Allocating z memory", FUNC_NAME, FAILURE);
    }

    ii = (int *) malloc (nums * sizeof (int));
    if (ii == NULL)
    {
        RETURN_ERROR ("Allocating ii memory", FUNC_NAME, FAILURE);
    }

    segments = nums - 1;
/* Songs original code
    if (nums % 4 != 0)
        segments++;
*/
    /* RDD - New code based on IDL */
    while (segments % 4 != 0)
        segments++;

    xmin = x[0];
    xmax = x[nums - 1];
    h = (xmax - xmin) / (float) segments;

    /* integrate spectral response over wavelength */
    /* Call spline to get second derivatives, */
    if (spline (x, f, nums, 2.0, 2.0, temp) != SUCCESS)
    {
        RETURN_ERROR ("Failed during spline", FUNC_NAME, FAILURE);
    }

    /* Call splint for interpolations. one-based arrays are considered */
    for (i = 0; i < nums; i++)
    {
        splint (x, f, temp, nums, i, &z[i]);
    }

    result = 0.0;
    /* Get the 5-points needed for Newton-Cotes formula */
#if 0
printf("-------------\n");
printf("nums = %d\n", nums);
#endif
    ii_count = (int) ((nums - 1) / 4);
    for (i = 0; i < ii_count; i++)
    {
        ii[i] = (i + 1) * 4;
#if 0
printf("ii[%d] = %d\n", i, ii[i]);
#endif
    }
    /* Compute the integral using the 5-point Newton-Cotes formula */
    for (i = 0; i < ii_count; i++)
    {
        result += 2.0 * h * (7.0 * (z[ii[i] - 4] + z[ii[i]]) +
                             32.0 * (z[ii[i] - 3] + z[ii[i] - 1]) +
                             12.0 * z[ii[i] - 2]) / 45.0;
    }

    *result_out = result;

    free (temp);
    free (z);
    free (ii);

    return SUCCESS;
}

/******************************************************************************
MODULE:  calculate_lt

PURPOSE: Calculate blackbody radiance from temperature using spectral repsonse 
         function

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/29/2014   Song Guo         Original Development
******************************************************************************/
int calculate_lt
(
    float temperature,         /*I: temperature */
    float **spectral_response, /*I: spectral response function */
    int num_srs,               /*I: number of spectral response points */
    float *radiance            /*O: blackbody radiance */
)
{
    char FUNC_NAME[] = "calculate_lt";
    int i;
    float rs_integral;
    float temp_integral;
    float *blackbody_radiance;
    float *product;

    /* Allocate memory */
    blackbody_radiance = (float *) malloc (num_srs * sizeof (float));
    if (blackbody_radiance == NULL)
    {
        RETURN_ERROR ("Allocating blackbody_radiance memory", FUNC_NAME,
                      FAILURE);
    }

    product = (float *) malloc (num_srs * sizeof (float));
    if (product == NULL)
    {
        RETURN_ERROR ("Allocating product memory", FUNC_NAME, FAILURE);
    }

    /* integrate spectral response over wavelength */
    if (int_tabulated (spectral_response[0], spectral_response[1], num_srs,
                       &rs_integral) != SUCCESS)
    {
        RETURN_ERROR ("Calling int_tabulated\n", FUNC_NAME, FAILURE);
    }

    /* using planck's equaiton to calculate radiance at each wavelength for 
       current temp */
    if (planck_eq (spectral_response[0], num_srs, temperature,
                   blackbody_radiance) != SUCCESS)
    {
        RETURN_ERROR ("Calling planck_eq\n", FUNC_NAME, FAILURE);
    }

    /* multiply the caluclated planck radiance by the spectral reponse and
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

    /* divide above result by integral of spectral response function */
    *radiance = temp_integral / rs_integral;

    /* Free allocated memory */
    free (blackbody_radiance);
    free (product);

    return SUCCESS;
}


/******************************************************************************
MODULE:  interpol

PURPOSE: Simulate IDL interpol function for lst purpose

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
10/3/2014   Song Guo         Original Development

TODO TODO TODO - RDD - I have modified this and have no idea (yet) if it is
                       correctly doing what is required of the IDL codebase.
                       But it was definitly wrong before I got it.
******************************************************************************/
void interpol
(
    float **v_x,     /*I: The input vector data */
    float *u,        /*I: The output grid points */
    int num_v,       /*I: number of input data */
    int num_u,       /*I: number of output grid points */
    int index,       /*I: Which index column data be used */
    float *r         /*O: Interpolated results */
)
{
    int i;
    int ix = 0;
    float s1;
    int d;

    if ((v_x[1][0] - v_x[0][0]) >= MINSIGMA)
    {
        s1 = 1;
    }
    else
    {
        s1 = -1;
    }

    for (i = 0; i < num_u; i++)
    {
        d = (int) (s1 * (u[i] - v_x[ix][0]));

        if (d == 0)
        {
            r[i] = v_x[ix][index];
        }
        else
        {
            if (d > 0)
            {
                while ((s1 * (u[i] - v_x[ix + 1][0]) > 0)
                       && (ix < num_v - 1))
                {
                    ix++;
                }
            }
            else
            {
                while ((s1 * (u[i] - v_x[ix][0]) < 0) && (ix > 0))
                {
                    ix--;
                }
            }

            r[i] = v_x[ix][index]
                   + (u[i] - v_x[ix][0]) * (v_x[ix + 1][index]
                                            - v_x[ix][index]) / (v_x[ix
                                                                 + 1][0]
                                                                 - v_x[ix][0]);
        }
    }
}


/******************************************************************************
MODULE:  calculate_lobs

PURPOSE: Calculate observed radiance from modtran results and the spectral
         reponse function

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/29/2014   Song Guo         Original Development
******************************************************************************/
int calculate_lobs
(
    float **modtran,           /*I: modtran results with wavelengths */
    float **spectral_response, /*I: spectral response function */
    int num_entries,           /*I: number of MODTRAN points */
    int num_srs,               /*I: number of spectral response points */
    int index,                 /*I: column index for data be used */
    float *radiance            /*O: LOB outputs */
)
{
    char FUNC_NAME[] = "calculate_lobs";
    int i;
    float *temp_rad;
    float rs_integral;
    float temp_integral;
    float *product;

    /* Allocate memory */
    temp_rad = (float *) malloc (num_entries * sizeof (float));
    if (temp_rad == NULL)
    {
        RETURN_ERROR ("Allocating temp_rad memory", FUNC_NAME, FAILURE);
    }

    product = (float *) malloc (num_entries * sizeof (float));
    if (product == NULL)
    {
        RETURN_ERROR ("Allocating product memory", FUNC_NAME, FAILURE);
    }

    /* integrate spectral response over wavelength */
    if (int_tabulated (spectral_response[0], spectral_response[1], num_srs,
                       &rs_integral) != SUCCESS)
    {
        RETURN_ERROR ("Calling int_tabulated\n", FUNC_NAME, FAILURE);
    }

    /* using planck's equaiton to calculate radiance at each wavelength for 
       current temp */
    interpol (modtran, spectral_response[0], num_entries, num_srs, index,
              temp_rad);

    /* multiply the caluclated radiance by the spectral reponse and integrate 
       over wavelength to get one number for current temp */
    for (i = 0; i < num_srs; i++)
    {
        product[i] = temp_rad[i] * spectral_response[1][i];
    }

    if (int_tabulated (spectral_response[0], product, num_srs,
                       &temp_integral) != SUCCESS)
    {
        RETURN_ERROR ("Calling int_tabulated\n", FUNC_NAME, FAILURE);
    }

    /* divide above result by integral of spectral response function */
    *radiance = temp_integral / rs_integral;

    /* Free allocated memory */
    free (temp_rad);
    free (product);

    return SUCCESS;
}


/******************************************************************************
MODULE:  calculate_point_atmospheric_parameters

PURPOSE: Generate transmission, upwelled radiance, and downwelled radiance at
         each height for each NARR point that is used.

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/29/2014   Song Guo         Original Development
******************************************************************************/
int calculate_point_atmospheric_parameters
(
    Input_t * input,       /*I: input structure */
    int num_points,        /*I: number of narr points */
    float albedo,          /*I: albedo */
    POINT_INFO *case_list, /*I: modtran run list */
    float **results,       /*O: atmospheric parameter for modtarn run */
    bool verbose           /*I: value to indicate if intermediate messages
                                should be printed */
)
{
    char FUNC_NAME[] = "calculate_point_atmospheric_parameters";
    FILE *fd;
    int i, j, k;
    int entry;
    float **spectral_response = NULL;
    float temp_radiance_0;
    float obs_radiance_0;
    float temp_radiance_273;
    float temp_radiance_300;
    int counter = 0;
    int index;
    int num_entries;
    int num_srs;       /* Number of spectral response values available */
    int result_loc;
    char current_file[PATH_MAX];
    float **temp1;
    float zero_temp;
    float **current_data;
    float x_0;
    float x_1;
    float inv_xx_diff; /* To save divisions */
    float y_0;
    float y_1;
    float tau; /* Transmission */
    float lu;  /* Upwelled Radiance */
    float ld;  /* Downwelled Radiance */
    /* TODO TODO TODO - This emissivity/albedo is just for water which is all
                        that has been implemented.  But the goal is to be
                        doing LAND, so integration aith the Aster data and
                        JPL conversion code will probably be needed before
                        execution of this routine, and then utilized here. */
    double emissivity = 1.0 - albedo;
    double inv_albedo = 1.0 / albedo;
    char *lst_data_dir = NULL;
    char full_path[PATH_MAX];
    int status;

    lst_data_dir = getenv ("LST_DATA_DIR");
    if (lst_data_dir == NULL)
    {
        RETURN_ERROR ("LST_DATA_DIR environment variable is not set",
                      FUNC_NAME, FAILURE);
    }

    LOG_MESSAGE ("Reading spectral response file", FUNC_NAME);
    if (input->meta.inst == INST_TM && input->meta.sat == SAT_LANDSAT_5)
    {
        num_srs = 171;

        /* Dynamic allocate the 2d memory */
        spectral_response =
            (float **) allocate_2d_array (2, num_srs, sizeof (float));
        if (spectral_response == NULL)
        {
            RETURN_ERROR ("Allocating spectral_response memory",
                          FUNC_NAME, FAILURE);
        }

        snprintf (full_path, sizeof (full_path),
                  "%s/%s", lst_data_dir, "L5_Spectral_Response.rsp");
        fd = fopen (full_path, "r");
        if (fd == NULL)
        {
            RETURN_ERROR ("Can't open L5_Spectral_Response.rsp file",
                          FUNC_NAME, FAILURE);
        }

        for (i = 0; i < num_srs; i++)
        {
            if (fscanf (fd, "%f %f%*c", &spectral_response[0][i],
                        &spectral_response[1][i]) != 2)
            {
                RETURN_ERROR ("Failed reading L5_Spectral_Response.rsp",
                              FUNC_NAME, FAILURE);
            }
        }
        fclose (fd);
    }
    else if (input->meta.inst == INST_ETM && input->meta.sat == SAT_LANDSAT_7)
    {
        num_srs = 47;

        /* Dynamic allocate the 2d memory */
        spectral_response =
            (float **) allocate_2d_array (2, num_srs, sizeof (float));
        if (spectral_response == NULL)
        {
            RETURN_ERROR ("Allocating spectral_response memory",
                          FUNC_NAME, FAILURE);
        }

        snprintf (full_path, sizeof (full_path),
                  "%s/%s", lst_data_dir, "L7_Spectral_Response.rsp");
        fd = fopen (full_path, "r");
        if (fd == NULL)
        {
            RETURN_ERROR ("Can't open L7_Spectral_Response.rsp file",
                          FUNC_NAME, FAILURE);
        }

        for (i = 0; i < num_srs; i++)
        {
            if (fscanf (fd, "%f %f%*c", &spectral_response[0][i],
                        &spectral_response[1][i]) != 2)
            {
                RETURN_ERROR ("Failed reading L7_Spectral_Response.rsp",
                              FUNC_NAME, FAILURE);
            }
        }
        fclose (fd);
    }
    else if (input->meta.inst == INST_OLI_TIRS
             && input->meta.sat == SAT_LANDSAT_8)
    {
        num_srs = 101;

        /* Dynamic allocate the 2d memory */
        spectral_response =
            (float **) allocate_2d_array (2, num_srs, sizeof (float));
        if (spectral_response == NULL)
        {
            RETURN_ERROR ("Allocating spectral_response memory",
                          FUNC_NAME, FAILURE);
        }

        snprintf (full_path, sizeof (full_path),
                  "%s/%s", lst_data_dir, "L8_B10.rsp");
        fd = fopen (full_path, "r");
        if (fd == NULL)
        {
            RETURN_ERROR ("Can't open L8_B10.rsp file", FUNC_NAME, FAILURE);
        }

        for (i = 0; i < num_srs; i++)
        {
            if (fscanf (fd, "%f %f%*c", &spectral_response[0][i],
                        &spectral_response[1][i]) == EOF)
            {
                RETURN_ERROR ("Failed reading L8_B10.rsp", FUNC_NAME, FAILURE);
            }
        }
        fclose (fd);
    }
    else
    {
        RETURN_ERROR ("invalid instrument type", FUNC_NAME, FAILURE);
    }

    /* calculate Lt for each temperature */
    if (input->meta.inst == INST_TM && input->meta.sat == SAT_LANDSAT_5)
    {
        if (calculate_lt (273, spectral_response, num_srs, &temp_radiance_273)
            != SUCCESS)
        {
            RETURN_ERROR ("Calling calculate_lt for 273K", FUNC_NAME, FAILURE);
        }
        if (calculate_lt (300, spectral_response, num_srs, &temp_radiance_300)
            != SUCCESS)
        {
            RETURN_ERROR ("Calling calculate_lt for 300K", FUNC_NAME, FAILURE);
        }
    }
    else if (input->meta.inst == INST_ETM && input->meta.sat == SAT_LANDSAT_7)
    {
        if (calculate_lt (273, spectral_response, num_srs, &temp_radiance_273)
            != SUCCESS)
        {
            RETURN_ERROR ("Calling calculate_lt for 273K", FUNC_NAME, FAILURE);
        }
        if (calculate_lt (300, spectral_response, num_srs, &temp_radiance_300)
            != SUCCESS)
        {
            RETURN_ERROR ("Calling calculate_lt for 300K", FUNC_NAME, FAILURE);
        }
    }
    else if (input->meta.inst == INST_OLI_TIRS
             && input->meta.sat == SAT_LANDSAT_8)
    {
        if (calculate_lt (273, spectral_response, num_srs, &temp_radiance_273)
            != SUCCESS)
        {
            RETURN_ERROR ("Calling calculate_lt for 273K", FUNC_NAME, FAILURE);
        }
        if (calculate_lt (300, spectral_response, num_srs, &temp_radiance_300)
            != SUCCESS)
        {
            RETURN_ERROR ("Calling calculate_lt for 300K", FUNC_NAME, FAILURE);
        }
    }
    else
    {
        RETURN_ERROR ("invalid instrument type", FUNC_NAME, FAILURE);
    }

    x_0 = temp_radiance_273;
    x_1 = temp_radiance_300;
    inv_xx_diff = 1.0F / (x_0 - x_1);

    /* iterate through all points and heights */
    for (i = 0; i < num_points; i++)
    {
        for (j = 0; j < NUM_ELEVATIONS; j++)
        {
            result_loc = i * NUM_ELEVATIONS + j;
            /* put results into results array */
            results[result_loc][LST_LATITUDE] = case_list[counter].latitude;
            results[result_loc][LST_LONGITUDE] = case_list[counter].longitude;
            results[result_loc][LST_HEIGHT] = case_list[counter].height;

            /* Read the lst_modtran.info file for the 000 execution
               We read the zero_temp from this file, and also the record count
               The record count is the same for all three associated runs */
            /* The 000 file is always the "counter+2" element in the array
               at this point in the code */
            snprintf (current_file, sizeof (current_file),
                      "%s/lst_modtran.info", case_list[counter+2].full_path);

            fd = fopen (current_file, "r");
            if (fd == NULL)
            {
                RETURN_ERROR ("Can't open current_file file",
                              FUNC_NAME, FAILURE);
            }
            /* determine temperature at lowest atmospheric layer
               (when MODTRAN is run at 0K) */
            if (fscanf (fd, "%*s %f%*c", &zero_temp) != 1)
            {
                RETURN_ERROR ("End of file (EOF) is met before"
                              " reading TARGET_PIXEL_SURFACE_TEMPERATURE",
                              FUNC_NAME, FAILURE);
            }
            /* determine number of entries in current file */
            if (fscanf (fd, "%*s %d%*c", &num_entries) != 1)
            {
                RETURN_ERROR ("End of file (EOF) is met before"
                              " reading RADIANCE_RECORD_COUNT",
                              FUNC_NAME, FAILURE);
            }
            fclose (fd);

            /* for each height, read in radiance inforomation for three
               modtran runs, columns of array are organized:
               wavelength | 273,0.0 | 310,0.0 | 000,0.1 */
            current_data =
                (float **) allocate_2d_array (num_entries, 4, sizeof (float));
            if (current_data == NULL)
            {
                RETURN_ERROR ("Allocating current_data memory",
                              FUNC_NAME, FAILURE);
            }

            temp1 =
                (float **) allocate_2d_array (num_entries, 2, sizeof (float));
            if (temp1 == NULL)
            {
                RETURN_ERROR ("Allocating temp1 memory", FUNC_NAME, FAILURE);
            }

            index = 0;
            /* iterate through three pairs of parameters */
            for (k = 0; k < 3; k++)
            {
                /* define current file */
                snprintf (current_file, sizeof (current_file),
                          "%s/lst_modtran.dat", case_list[counter].full_path);

                fd = fopen (current_file, "r");
                if (fd == NULL)
                {
                    RETURN_ERROR ("Can't open current_file file",
                                  FUNC_NAME, FAILURE);
                }
                for (entry = 0; entry < num_entries; entry++)
                {
                    if (fscanf (fd, "%f %f%*c", &temp1[entry][0],
                                &temp1[entry][1]) != 2)
                    {
                        RETURN_ERROR ("Failed reading lst_modtran.dat lines",
                                      FUNC_NAME, FAILURE);
                    }
                }
                fclose (fd);

                /* put arrays into data array for current point at current
                   height */
                if (index == 0)
                {
                    for (entry = 0; entry < num_entries; entry++)
                    {
                        current_data[entry][0] = temp1[entry][0];
                        current_data[entry][1] = temp1[entry][1];
                    }
                    index++;
                }
                else
                {
                    for (entry = 0; entry < num_entries; entry++)
                    {
                        current_data[entry][index] = temp1[entry][1];
                    }
                }

                counter++;
                index++;
            }

#if 0
            printf ("num_srs = %d\n", num_srs);
            printf ("num_entries = %d\n", num_entries);
            printf ("temp_radiance_273 = %f\n", temp_radiance_273);
            printf ("temp_radiance_300 = %f\n", temp_radiance_300);
            fflush(stdout);
#endif
            /* parameters from 3 modtran runs
               Lobs = Lt*tau + Lu; m = tau; b = Lu; */
            status = calculate_lobs (current_data, spectral_response,
                                     num_entries, num_srs, 1, &y_0);
            if (status != SUCCESS)
            {
                RETURN_ERROR ("Calling calculate_lob 1", FUNC_NAME, FAILURE);
            }

            status = calculate_lobs (current_data, spectral_response,
                                     num_entries, num_srs, 2, &y_1);
            if (status != SUCCESS)
            {
                RETURN_ERROR ("Calling calculate_lob 2", FUNC_NAME, FAILURE);
            }

            if (free_2d_array ((void **) temp1) != SUCCESS)
            {
                RETURN_ERROR ("Freeing memory: temp\n", FUNC_NAME, FAILURE);
            }

            /* Implement a = INVERT(TRANSPOSE(x)##x)##TRANSPOSE(x)##y 
               Note: I slove the two equations analytically */
            tau = (y_0 - y_1) * inv_xx_diff;
            lu = (x_1 * y_0 - x_0 * y_1) * inv_xx_diff;

            /* determine Lobs and Lt when
               modtran was run at 0K - calculate downwelled */
            if (calculate_lt (zero_temp, spectral_response, num_srs,
                              &temp_radiance_0) != SUCCESS)
            {
                RETURN_ERROR ("Calling calculate_lt for 0K",
                              FUNC_NAME, FAILURE);
            }

            if (calculate_lobs (current_data, spectral_response,
                                num_entries, num_srs, 3, &obs_radiance_0)
                != SUCCESS)
            {
                RETURN_ERROR ("Calling calculate_lob 2", FUNC_NAME, FAILURE);
            }

            /* Free the allocated memory in the loop */
            if (free_2d_array ((void **) current_data) != SUCCESS)
            {
                RETURN_ERROR ("Freeing memory: current_data\n",
                              FUNC_NAME, FAILURE);
            }

            /* Ld = (((Lobs-Lu)/tau) - (Lt*emissivity))/(1.0-emissivity) */
            /* Ld = (((Lobs-Lu)/tau) - (Lt*emissivity))/abledo */
            /* Ld = (((Lobs-Lu)/tau) - (Lt*emissivity))*inv_albedo */
            ld = (((obs_radiance_0 - lu) / tau)
                  - (temp_radiance_0 * emissivity)) * inv_albedo;

            /* put remaining results into results array */
            results[result_loc][LST_TRANSMISSION] = tau;
            results[result_loc][LST_UPWELLED_RADIANCE] = lu;
            results[result_loc][LST_DOWNWELLED_RADIANCE] = ld;
        } /* END - NUM_ELEVATIONS loop */
    } /* END - num_points loop */

    /* Free allocated memory */
    if (free_2d_array ((void **) spectral_response) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: spectral_response\n", FUNC_NAME,
                      FAILURE);
    }

    /* write results to a file */
    snprintf (current_file, sizeof (current_file),
              "atmosphericParameters.txt");
    printf ("Creating Atmospheric Parameters File = [%s]\n", current_file);
    fd = fopen (current_file, "w");
    if (fd == NULL)
    {
        RETURN_ERROR ("Can't open atmosphericParameters.txt file",
                      FUNC_NAME, FAILURE);
    }
    for (k = 0; k < num_points * NUM_ELEVATIONS; k++)
    {
        fprintf (fd, "%f,%f,%f,%f,%f,%f\n",
                 results[k][LST_LATITUDE],
                 results[k][LST_LONGITUDE],
                 results[k][LST_HEIGHT],
                 results[k][LST_TRANSMISSION],
                 results[k][LST_UPWELLED_RADIANCE],
                 results[k][LST_DOWNWELLED_RADIANCE]);
    }
    fclose (fd);

    return SUCCESS;
}
