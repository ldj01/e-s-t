/*****************************************************************************
NAME: integrate.c

PURPOSE: This file contains a set of functions associated with tabular
         integration.
*****************************************************************************/
#include <stdlib.h>
#include "const.h"
#include "utilities.h"
#include "integrate.h"

static integration_workspace workspace = {0, {NULL, NULL}, 0, NULL, 0, NULL};

static double **get_spline_workspace(int size)
{
    if (workspace.size < size)
    {
        workspace.array[0] = realloc(workspace.array[0], size*sizeof(double));
        workspace.array[1] = realloc(workspace.array[1], size*sizeof(double));
        if (workspace.array[0] == NULL || workspace.array[1] == NULL)
            free_integration_workspace();
        else
            workspace.size = size;
    }

    return workspace.array;
}

static double *get_integration_z(int size)
{
    if (workspace.zsize < size)
    {
        workspace.z = realloc(workspace.z, size*sizeof(double));
        if (workspace.z == NULL)
            free_integration_workspace();
        else
            workspace.zsize = size;
    }

    return workspace.z;
}

static int *get_integration_i(int size)
{
    if (workspace.isize < size)
    {
        workspace.i = realloc(workspace.i, size*sizeof(int));
        if (workspace.i == NULL)
            free_integration_workspace();
        else
            workspace.isize = size;
    }

    return workspace.i;
}

void free_integration_workspace()
{
    free(workspace.array[0]);
    free(workspace.array[1]);
    workspace.array[0] = NULL;
    workspace.array[1] = NULL;
    workspace.size = 0;

    free(workspace.z);
    workspace.z = NULL;
    workspace.zsize = 0;

    free(workspace.i);
    workspace.i = NULL;
    workspace.isize = 0;
}


/*****************************************************************************
MODULE:  spline

PURPOSE: Natural spline constructs a cubic spline through a set of x and y
         values.

RETURN: SUCCESS
        FAILURE
*****************************************************************************/
static int spline
(
    double *x,
    double *y,
    int n,
    double *y2
)
{
    char FUNC_NAME[] = "spline";
    int i;
    double p;
    double qn;
    double sig;
    double un;
    double *u = get_spline_workspace(n)[0];

    if (u == NULL)
    {
        RETURN_ERROR ("Can't allocate memory", FUNC_NAME, FAILURE);
    }

    /* Set the lower boundary */
    y2[0] = 0.0;
    u[0] = 0.0;

    /* Set the upper boundary */
    qn = 0.0;
    un = 0.0;

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

    return SUCCESS;
}


/*****************************************************************************
MODULE:  splint

PURPOSE: Use the cubic spline generated with spline to interpolate
         values in the XY table.
*****************************************************************************/
static void splint
(
    double *xa,
    double *ya,
    double *y2a,
    int n,
    double x,
    double *y
)
{
    int k;
    double h;
    double a;
    static int splint_klo = -1;
    static int splint_khi = -1;
    static double one_sixth = (1.0 / 6.0); /* To remove a division */

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

        /*  The equation used below is the following, simplified:

            b = 1 - a;
            *y = a * ya[splint_klo]
               + b * ya[splint_khi]
               + ((a * a * a - a) * y2a[splint_klo]
               + (b * b * b - b) * y2a[splint_khi]) * (h * h) * one_sixth; */

        *y = ya[splint_khi] + a*(ya[splint_klo] - ya[splint_khi])
            + one_sixth*h*h*a*(a - 1)*((a + 1)*y2a[splint_klo] +
                                       (2 - a)*y2a[splint_khi]);
    }
}


/*****************************************************************************
MODULE:  int_tabulated

PURPOSE: This function integrates a tabulated set of data { x(i) , f(i) },
         on the closed interval [min(X) , max(X)].

RETURN: SUCCESS
        FAILURE

NOTE: x and f are assumed to be in sorted order (min(x) -> max(x))
*****************************************************************************/
int int_tabulated
(
    double *x,         /*I: Tabulated X-value data */
    double *f,         /*I: Tabulated F-value data */
    int nums,          /*I: Number of points */
    double *result_out /*O: Integrated result */
)
{
    char FUNC_NAME[] = "int_tabulated";
    double *temp = NULL;
    double *z = NULL;
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

    /* Get the working space arrays. */
    temp = get_spline_workspace(nums)[1];
    if (temp == NULL)
    {
        RETURN_ERROR ("Allocating temp memory", FUNC_NAME, FAILURE);
    }

    z = get_integration_z(segments + 1);
    if (z == NULL)
    {
        RETURN_ERROR ("Allocating z memory", FUNC_NAME, FAILURE);
    }

    ii = get_integration_i(ii_count);
    if (ii == NULL)
    {
        RETURN_ERROR ("Allocating ii memory", FUNC_NAME, FAILURE);
    }

    /* Interpolate spectral response over wavelength */
    if (spline (x, f, nums, temp) != SUCCESS)
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
        double *z_ptr = &z[ii[i] - 4];

        result += 14*(z_ptr[0] + z_ptr[4]) + 64*(z_ptr[1] + z_ptr[3])
                + 24*z_ptr[2];
    }

    /* Assign the results to the output */
    *result_out = result*h/45;

    return SUCCESS;
}
