#include <stdarg.h>
#include <math.h>


#include "const.h"
#include "2d_array.h"
#include "utilities.h"
#include "input.h"
#include "lst_types.h"


/******************************************************************************
METHOD:  planck_eq

PURPOSE: Using Planck's equaiton to calculate radiance at each wavelength for
         current temperature.

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/30/2014   Song Guo         Original Development
******************************************************************************/
void planck_eq
(
    float *wavelength,  /* I: Each wavelength */
    int num_elements,   /* I: Number of wavelengths to calculate */
    float temperature,  /* I: The temperature to calculate for */
    double *bb_radiance /* I: the blackbody results for each wavelength */
)
{
    int i;
    double lambda;

    /* Planck Const hecht pg, 585 ## units: Js */
    double PLANCK_CONST = (6.6260755 * pow (10, -34));

    /* Boltzmann Gas Const halliday et 2001 -- units: J/K */
    double BOLTZMANN_GAS_CONST = (1.3806503 * pow (10, -23));

    /* Speed of Light -- units: m/s */
    double SPEED_OF_LIGHT = (299792458.0);
    double SPEED_OF_LIGHT_SQRD = (SPEED_OF_LIGHT * SPEED_OF_LIGHT);

    for (i = 0; i < num_elements; i++)
    {
        /* Lambda intervals of spectral response locations microns units: m */
        lambda = wavelength[i] * pow (10, -6);

        /* Compute the Planck Blackbody Eq [W/m^2 sr um] */
        bb_radiance[i] = 2.0 * PLANCK_CONST * SPEED_OF_LIGHT_SQRD
                         * (pow (10, -6) * pow (lambda, -5.0))
                         * (1.0 / (exp ((PLANCK_CONST * SPEED_OF_LIGHT)
                                         / (lambda
                                            * BOLTZMANN_GAS_CONST
                                            * temperature))
                                   - 1.0));

        /* convert to W/cm^2 sr micron to match modtran units */
        /* br / (100 * 100) == br * 10e-5 */
        bb_radiance[i] *= 10e-5;
    }
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

    u = malloc ((unsigned) (n - 1) * sizeof (double));
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
static int splint_klo = -1;
static int splint_khi = -1;
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

    if (splint_klo < 0)
    {
        splint_klo = 0;
        splint_khi = n - 1;
    }
    else
    {
        if (x < xa[splint_klo])
            splint_klo = 0;
        if (x > xa[splint_khi])
            splint_khi = n - 1;
    }

    while (splint_khi - splint_klo > 1)
    {
        k = (splint_khi + splint_klo) >> 1;

        if (xa[k] > x)
            splint_khi = k;
        else
            splint_klo = k;
    }

    h = xa[splint_khi] - xa[splint_klo];

    if (h == 0.0)
    {
        *y = 0.0;
    }
    else
    {
        a = (xa[splint_khi] - x) / h;

        b = (x - xa[splint_klo]) / h;

        *y = a * ya[splint_klo]
             + b * ya[splint_khi]
             + ((a * a * a - a) * y2a[splint_klo]
                + (b * b * b - b) * y2a[splint_khi]) * (h * h) * one_sixth;
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

NOTE: x and f are assumed to be in sorted order (min(x) -> max(x))
******************************************************************************/
int int_tabulated
(
    float *x,         /*I: Tabulated X-value data */
    float *f,         /*I: Tabulated F-value data */
    int nums,         /*I: Number of points */
    float *result_out /*O: Integrated result */
)
{
    char FUNC_NAME[] = "int_tabulated";
    float *temp = NULL;
    float *z = NULL;
    double xmin;
    double xmax;
    int i;
    int *ii = NULL;
    int ii_count;
    double h;
    double result;
    int segments;

    /* Figure out the number of segments needed */
    segments = nums - 1;
    while (segments % 4 != 0)
        segments++;

    /* Determine how many iterations are needed  */
    ii_count = (int) ((segments) / 4);

    /* Determine the min and max */
    xmin = x[0];
    xmax = x[nums - 1];

    /* Determine the step size */
    h = (xmax - xmin) / segments;

    /* Allocate memory */
    temp = malloc (nums * sizeof (float));
    if (temp == NULL)
    {
        RETURN_ERROR ("Allocating temp memory", FUNC_NAME, FAILURE);
    }

    z = malloc ((segments+1) * sizeof (float));
    if (z == NULL)
    {
        RETURN_ERROR ("Allocating z memory", FUNC_NAME, FAILURE);
    }

    ii = malloc (ii_count * sizeof (int));
    if (ii == NULL)
    {
        RETURN_ERROR ("Allocating ii memory", FUNC_NAME, FAILURE);
    }

    /* Interpolate spectral response over wavelength */
    /* Using 1e30 forces generation of a natural spline and produces nearly
       the same results as IDL */
    if (spline (x, f, nums, 1e30, 1e30, temp) != SUCCESS)
    {
        RETURN_ERROR ("Failed during spline", FUNC_NAME, FAILURE);
    }

    /* Call splint for interpolations. one-based arrays are considered */
    for (i = 0; i < segments+1; i++)
    {
        splint (x, f, temp, nums, h*i+xmin, &z[i]);
    }

    /* Get the 5-points needed for Newton-Cotes formula */
    for (i = 0; i < ii_count; i++)
    {
        ii[i] = (i + 1) * 4;
    }

    /* Compute the integral using the 5-point Newton-Cotes formula */
    result = 0.0;
    for (i = 0; i < ii_count; i++)
    {
        result += (h * (14.0 * (z[ii[i] - 4] + z[ii[i]]) +
                        64.0 * (z[ii[i] - 3] + z[ii[i] - 1]) +
                        24.0 * z[ii[i] - 2]) / 45.0);
    }

    /* Assign the results to the output */
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
    float *product;
    double *blackbody_radiance;

    /* Allocate memory */
    blackbody_radiance = malloc (num_srs * sizeof (double));
    if (blackbody_radiance == NULL)
    {
        RETURN_ERROR ("Allocating blackbody_radiance memory", FUNC_NAME,
                      FAILURE);
    }

    product = malloc (num_srs * sizeof (float));
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

    /* Use planck's blackbody radiance equation to calculate radiance at each
       wavelength for the current temperature */
    planck_eq (spectral_response[0], num_srs, temperature, blackbody_radiance);

    /* multiply the calculated planck radiance by the spectral reponse and
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
MODULE:  linear_interpolate_over_modtran

PURPOSE: Simulate IDL (interpol) function for LST.

RETURN: SUCCESS
        FAILURE
******************************************************************************/
void linear_interpolate_over_modtran
(
    float **modtran, /* I: The MODTRAN data - provides both the a and b */
    int index,       /* I: The MODTRAN temperatur to use for a */
    float *c,        /* I: The Landsat wavelength grid points */
    int num_in,      /* I: Number of input data and grid points*/
    int num_out,     /* I: Number of output grid points */
    float *x         /* O: Interpolated output results */
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
        x[o] = d1 + ((g - g1) / (g2 - g1)) * (d2 - d1);
    }
}


/******************************************************************************
MODULE:  calculate_lobs

PURPOSE: Calculate observed radiance from MODTRAN results and the spectral
         reponse function.

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/29/2014   Song Guo         Original Development
******************************************************************************/
int calculate_lobs
(
    float **modtran,           /*I: MODTRAN results with wavelengths */
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
    temp_rad = malloc (num_srs * sizeof (float));
    if (temp_rad == NULL)
    {
        RETURN_ERROR ("Allocating temp_rad memory", FUNC_NAME, FAILURE);
    }

    product = malloc (num_srs * sizeof (float));
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

    /* interpolate MODTRAN radiance to Landsat wavelengths */
    linear_interpolate_over_modtran (modtran, index, spectral_response[0],
                                     num_entries, num_srs, temp_rad);

    /* multiply the calculated radiance by the spectral reponse and integrate
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

    /* divide above result by integral of spectral response function */
    *radiance = temp_integral / rs_integral;

    /* Free allocated memory */
    free (temp_rad);
    free (product);

    return SUCCESS;
}


/*****************************************************************************
METHOD:  matrix_transpose_2x2

PURPOSE: Transposes a 2x2 matrix, producing a 2x2 result.

RETURN: None
*****************************************************************************/
void matrix_transpose_2x2(double *A, double *out)
{
    /*
        Formula is:

            out[0] = a
            out[1] = c
            out[2] = b
            out[3] = d

        Where:

            a = A[0]
            b = A[1]
            c = A[2]
            d = A[3]
    */

    out[0] = A[0];
    out[1] = A[2];
    out[2] = A[1];
    out[3] = A[3];
}


/*****************************************************************************
METHOD:  matrix_inverse_2x2

PURPOSE: Inverts a 2x2 matrix, producing a 2x2 result.

RETURN: None
*****************************************************************************/
void matrix_inverse_2x2(double *A, double *out)
{
    /*
        Formula is:

            out[0] = d * determinant;
            out[1] = (-b) * determinant;
            out[2] = (-c) * determinant;
            out[3] = a * determinant;

        Where:

            a = A[0]
            b = A[1]
            c = A[2]
            d = A[3]

            determinant = 1.0 / (a * d - b * c)
    */

    double determinant = (1.0 / (A[0] * A[3] - A[1] * A[2]));

    out[0] = A[3] * determinant;
    out[1] = (-A[1]) * determinant;
    out[2] = (-A[2]) * determinant;
    out[3] = A[0] * determinant;
}


/*****************************************************************************
METHOD:  matrix_multiply_2x2_2x2

PURPOSE: Multiply a 2x2 matrix with a 2x2 matrix, producing a 2x2 result.

RETURN: None
*****************************************************************************/
void matrix_multiply_2x2_2x2(double *A, double *B, double *out)
{
    /*
        Formula is:

            out[0] = a * e + b * g
            out[1] = a * f + b * h
            out[2] = c * e + d * g
            out[3] = c * f + d * h

        Where:

            a = A[0]
            b = A[1]
            c = A[2]
            d = A[3]

            e = B[0]
            f = B[1]
            g = B[2]
            h = B[3]
    */

    out[0] = A[0] * B[0] + A[1] * B[2];
    out[1] = A[0] * B[1] + A[1] * B[3];
    out[2] = A[2] * B[0] + A[3] * B[2];
    out[3] = A[2] * B[1] + A[3] * B[3];
}


/*****************************************************************************
METHOD:  matrix_multiply_2x2_2x1

PURPOSE: Multiply a 2x2 matrix with a 2x1 matrix, producing a 2x1 result.

RETURN: None
*****************************************************************************/
void matrix_multiply_2x2_2x1(double *A, double *B, double *out)
{
    /*
        Formula is:

            out[0] = a * e + b * f
            out[1] = c * e + d * f

        Where:

            a = A[0]
            b = A[1]
            c = A[2]
            d = A[3]

            e = B[0]
            f = B[1]
    */

    out[0] = A[0] * B[0] + A[1] * B[1];
    out[1] = A[2] * B[0] + A[3] * B[1];
}


/*****************************************************************************
METHOD:  calculate_point_atmospheric_parameters

PURPOSE: Generate transmission, upwelled radiance, and downwelled radiance at
         each height for each NARR point that is used.

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/29/2014   Song Guo         Original Development
*****************************************************************************/
#define L5_TM_SRS_COUNT (171)
#define L7_TM_SRS_COUNT (47)
#define L8_OLITIRS_SRS_COUNT (101)
#define MAX_SRS_COUNT (L5_TM_SRS_COUNT)
int calculate_point_atmospheric_parameters
(
    Input_t *input,            /* I: input structure */
    REANALYSIS_POINTS *points, /* I: The coordinate points */
    float albedo,              /* I: albedo */
    float **modtran_results,   /* O: atmospheric parameter for modtarn run */
    bool verbose               /* I: value to indicate if intermediate
                                     messages should be printed */
)
{
    char FUNC_NAME[] = "calculate_point_atmospheric_parameters";

    FILE *fd;
    FILE *used_points_fd;

    int i;
    int j;
    int k;
    int entry;

    float **spectral_response = NULL;
    float temp_radiance_0;
    float obs_radiance_0;
    float temp_radiance_273;
    float temp_radiance_310;
    int counter;
    int index;
    int num_entries;   /* Number of MODTRAN output results to read and use */
    int num_srs;       /* Number of spectral response values available */
    int result_loc;

    char *lst_data_dir = NULL;
    char current_file[PATH_MAX];
    char srs_file_path[PATH_MAX];
    char msg[PATH_MAX];

    float modtran_wavelength;
    float modtran_radiance;
    float zero_temp;
    float **current_data;
    float y_0;
    float y_1;
    float tau; /* Transmission */
    float lu;  /* Upwelled Radiance */
    float ld;  /* Downwelled Radiance */
    /* TODO TODO TODO - This emissivity/albedo is just for water which is all
                        that has been implemented.  But the goal is to be
                        doing LAND, so integration with the Aster data and
                        JPL conversion code will probably be needed before
                        execution of this routine, and then utilized here. */
    double emissivity = 1.0 - albedo;
    double inv_albedo = 1.0 / albedo;

    /* Variables to hold matricies and the results for the operations perfomed
       on them */
    double X_2x2[4];
    double Xt_2x2[4];
    double Xt_X_2x2[4];
    double Inv_Xt_X_2x2[4];
    double Y_2x1[2];
    double Xt_Y_2x1[4];
    double A_2x1[2];


    lst_data_dir = getenv ("LST_DATA_DIR");
    if (lst_data_dir == NULL)
    {
        RETURN_ERROR ("LST_DATA_DIR environment variable is not set",
                      FUNC_NAME, FAILURE);
    }

    /* Allocate memory for maximum spectral response count */
    spectral_response =
        (float **) allocate_2d_array (2, MAX_SRS_COUNT, sizeof (float));
    if (spectral_response == NULL)
    {
        RETURN_ERROR ("Allocating spectral_response memory",
                      FUNC_NAME, FAILURE);
    }

    /* Determine the spectral response file to read */
    if (input->meta.instrument == INST_TM
        && input->meta.satellite == SAT_LANDSAT_5)
    {
        num_srs = L5_TM_SRS_COUNT;

        snprintf (srs_file_path, sizeof (srs_file_path),
                  "%s/%s", lst_data_dir, "L5_Spectral_Response.rsp");
    }
    else if (input->meta.instrument == INST_ETM
             && input->meta.satellite == SAT_LANDSAT_7)
    {
        num_srs = L7_TM_SRS_COUNT;

        snprintf (srs_file_path, sizeof (srs_file_path),
                  "%s/%s", lst_data_dir, "L7_Spectral_Response.rsp");
    }
    else if (input->meta.instrument == INST_OLI_TIRS
             && input->meta.satellite == SAT_LANDSAT_8)
    {
        num_srs = L8_OLITIRS_SRS_COUNT;

        snprintf (srs_file_path, sizeof (srs_file_path),
                  "%s/%s", lst_data_dir, "L8_B10.rsp");
    }
    else
    {
        RETURN_ERROR ("invalid instrument type", FUNC_NAME, FAILURE);
    }

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
        if (fscanf (fd, "%f %f%*c", &spectral_response[0][i],
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

    /* Implement a = INVERT(TRANSPOSE(x)##x)##TRANSPOSE(x)##y
       from the IDL code base.
       Partially implemented here, the variable part is implemented in the
       looping code. */
    X_2x2[0] = 1;
    X_2x2[1] = temp_radiance_273;
    X_2x2[2] = 1;
    X_2x2[3] = temp_radiance_310;

    matrix_transpose_2x2(X_2x2, Xt_2x2);
    matrix_multiply_2x2_2x2(Xt_2x2, X_2x2, Xt_X_2x2);
    matrix_inverse_2x2(Xt_X_2x2, Inv_Xt_X_2x2);

    /* Output information about the used points, primarily usefull for
       plotting them against the scene */
    used_points_fd = fopen ("used_points.txt", "w");
    if (used_points_fd == NULL)
    {
        RETURN_ERROR ("Can't open used_points.txt file",
                      FUNC_NAME, FAILURE);
    }

    /* Iterate through all points and heights */
    counter = 0;
    for (i = 0; i < points->num_points; i++)
    {
        fprintf (used_points_fd, "\"%d\"|\"%f\"|\"%f\"\n",
                 i, points->utm_easting[i], points->utm_northing[i]);

        for (j = 0; j < NUM_ELEVATIONS; j++)
        {
            result_loc = i * NUM_ELEVATIONS + j;

            /* put results into MODTRAN results array */
            modtran_results[result_loc][LST_LATITUDE] =
                points->modtran_runs[counter].latitude;
            modtran_results[result_loc][LST_LONGITUDE] =
                points->modtran_runs[counter].longitude;
            modtran_results[result_loc][LST_HEIGHT] =
                points->modtran_runs[counter].height;

            /* Read the lst_modtran.info file for the 000 execution
               (when MODTRAN is run at 0K)
               We read the zero_temp from this file, and also the record count
               The record count is the same for all three associated runs */
            /* The 000 file is always the "counter+2" element in the array
               at this point in the code */
            snprintf (current_file, sizeof (current_file),
                      "%s/lst_modtran.info",
                      points->modtran_runs[counter+2].path);

            fd = fopen (current_file, "r");
            if (fd == NULL)
            {
                RETURN_ERROR ("Can't open current_file file",
                              FUNC_NAME, FAILURE);
            }
            /* Retrieve the temperature from this lowest atmospheric layer */
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

            /* iterate through the three pairs of parameters */
            for (index = 1; index < 4; index++)
            {
                /* define current file */
                snprintf (current_file, sizeof (current_file),
                          "%s/lst_modtran.dat",
                          points->modtran_runs[counter].path);

                fd = fopen (current_file, "r");
                if (fd == NULL)
                {
                    RETURN_ERROR ("Can't open current_file file",
                                  FUNC_NAME, FAILURE);
                }
                for (entry = 0; entry < num_entries; entry++)
                {
                    if (fscanf (fd, "%f %f%*c",
                                &modtran_wavelength, &modtran_radiance)
                        != 2)
                    {
                        RETURN_ERROR ("Failed reading lst_modtran.dat lines",
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

            /* parameters from 3 modtran runs
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

            /* Implement a = INVERT(TRANSPOSE(x)##x)##TRANSPOSE(x)##y
               from the IDL code base.
               Partially implemented above and used here. */
            Y_2x1[0] = y_0;
            Y_2x1[1] = y_1;

            matrix_multiply_2x2_2x1(Xt_2x2, Y_2x1, Xt_Y_2x1);
            matrix_multiply_2x2_2x1(Inv_Xt_X_2x2, Xt_Y_2x1, A_2x1);

            tau = A_2x1[1]; /* Transmittance */
            lu = A_2x1[0];  /* Upwelled Radiance */

            /* determine Lobs and Lt when
               modtran was run at 0K - calculate downwelled */
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

            /* Calculate the downwelled radiance */
            /* TODO - These are all equivalent today, but may need to change
                      to use the (1.0 - emissivity) when emissivity is
                      provided. */
            /* Ld = (((Lobs-Lu)/tau) - (Lt*emissivity))/(1.0-emissivity) */
            /* Ld = (((Lobs-Lu)/tau) - (Lt*emissivity))/abledo */
            /* Ld = (((Lobs-Lu)/tau) - (Lt*emissivity))*inv_albedo */
            ld = (((obs_radiance_0 - lu) / tau)
                  - (temp_radiance_0 * emissivity)) * inv_albedo;

            /* Place results into MODTRAN results array */
            modtran_results[result_loc][LST_TRANSMISSION] = tau;
            modtran_results[result_loc][LST_UPWELLED_RADIANCE] = lu;
            modtran_results[result_loc][LST_DOWNWELLED_RADIANCE] = ld;

            /* Free the allocated memory in the loop */
            if (free_2d_array ((void **) current_data) != SUCCESS)
            {
                RETURN_ERROR ("Freeing memory: current_data\n",
                              FUNC_NAME, FAILURE);
            }
            current_data = NULL;
        } /* END - NUM_ELEVATIONS loop */
    } /* END - num_points loop */
    fclose (used_points_fd);

    /* Free allocated memory */
    if (free_2d_array ((void **) spectral_response) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: spectral_response\n", FUNC_NAME,
                      FAILURE);
    }
    spectral_response = NULL;

    /* Output the results to a file */
    snprintf (current_file, sizeof (current_file),
              "atmosphericParameters.txt");
    snprintf (msg, sizeof (msg),
              "Creating Atmospheric Parameters File = [%s]\n", current_file);
    LOG_MESSAGE (msg, FUNC_NAME);
    fd = fopen (current_file, "w");
    if (fd == NULL)
    {
        RETURN_ERROR ("Can't open atmosphericParameters.txt file",
                      FUNC_NAME, FAILURE);
    }
    for (k = 0; k < points->num_points * NUM_ELEVATIONS; k++)
    {
        fprintf (fd, "%f,%f,%12.9f,%12.9f,%12.9f,%12.9f\n",
                 modtran_results[k][LST_LATITUDE],
                 modtran_results[k][LST_LONGITUDE],
                 modtran_results[k][LST_HEIGHT],
                 modtran_results[k][LST_TRANSMISSION],
                 modtran_results[k][LST_UPWELLED_RADIANCE],
                 modtran_results[k][LST_DOWNWELLED_RADIANCE]);
    }
    fclose (fd);

    return SUCCESS;
}
