#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "2d_array.h"
#include "input.h"
#include "const.h"
#include "scene_based_lst.h"

static int klo = -1 ;
static int khi = -1 ;

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
void planck_eq
(
    float *wavelength,
    int num_srs, 
    float temperature,
    float *black_radiance
)
{
    int i;
    float *lambda;
    float h;
    float k;
    float c;
   
    /* Allocate memory */
    lambda = (float *)malloc(num_srs * sizeof(float));  
    if (temp == NULL)
    {
        sprintf (errstr, "Allocating lambda memory");
        LST_ERROR (errstr, "second_narr");
    }

    /* Lambda intervals of Landsat5 spectral response locations microns 
       units: m */
    lambda = wavelengths * (float)10e-6;

    /* Planck Const hecht pg, 585 ## units: Js */
    h = 6.6260755 * (float)10e-34;

    /* Boltzmann Gas Const halliday et 2001 ## units: J/K */
    k = 1.3806503 * (float)10e-23;

    /* Speed of Light ## units: m/s */
    c = 299792458.0

     for (i = 0; i < num_srs; i++)
    {
        /* Lambda intervals of Landsat5 spectral response locations microns 
           units: m */
        lambda[i] = wavelengths[i] * (float)10e-6;

        /* Compute the Planck Blackbody Eq [W/m^2 sr um] */  
        blank_radiance[i] = 2.0 * h * (pow(c, 2.0)) * ((float)(10e-6) * 
            pow(lambda[i],-5.0)) * (1.0/(exp((h * c) / (lambda[i] * k * 
            temperature)) -1.0));

        /* convert to W/cm^2 sr micron to match modtran units */
        blank_radiance[i] /= 10000.0; 
    }
    free(lambda);

    return SUCCESS;
}

/******************************************************************************
MODULE:  calculate_lobs

PURPOSE: Calculate observed radiance from modtran results and the spectral reponse 
         function

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/29/2014   Song Guo         
******************************************************************************/
int calculate_lobs
(
    float wavelengths,          /*I: wavelength */
    float *modtran,             /*I: modtran results */
    float **spectral_response,  /*I: spectral response function */
    int num_srs,                /*I: number of spectral response points */
    float *radiance             /*O: LOB outputs */  
)
{ 
    int i;
    float *temp_rad;
    float rsr_integral;
    float temp_integral;
    float *product;

    /* Allocate memory */
    temp_rad = (float *)malloc(num_srs * sizeof(float));  
    if (temp_rad == NULL)
    {
        sprintf (errstr, "Allocating temp_rad memory");
        LST_ERROR (errstr, "second_narr");
    }

    product = (float *)malloc(num_srs * sizeof(float));  
    if (product == NULL)
    {
        sprintf (errstr, "Allocating product memory");
        LST_ERROR (errstr, "second_narr");
    }

    /* integrate spectral response over wavelength */
    int_tabulated(spectral_response[0], spectral_response[1], &rs_integral);

    /* using planck's equaiton to calculate radiance at each wavelength for 
       current temp*/
    status = interpol(modtran, wavelengths, spectral_response[0], num_srs, 
                      temp_rad);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Calling interpol\n");
        RETURN_ERROR (errstr, "second_narr", FAILURE);              
    }

    /* multiply the caluclated radiance by the spectral reponse and integrate 
       over wavelength to get one number for current temp*/
    for (i = 0; i < num_srs; i++) 
    {
        product[i] = temp_rad[i] * spectral_response[1][i];
    }

    int_tabulated(spectral_response[0], product, &temp_integral);

    /* divide above result by integral of spectral response function */
    radiance = temp_integral / rs_integral;

    /* Free allocated memory */
    free(temp_rad);
    free(product);

    return SUCCESS;

}

/******************************************************************************
MODULE:  third_pixels_post

PURPOSE: Generate transmission, upwelled radiance, and downwelled radiance at 
         each Landsat pixel

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
10/7/2014   Song Guo         Original Development
******************************************************************************/
int third_pixels_post
(
    Input_t *input,             /*I: input structure */
    int num_points,             /*I: number of narr points */
    char **dem_infile,          /*I: address of input DEM filename */
    char **emi_infile,          /*I: address of input ASTER emissivity filename */
    double **results,           /*I: atmospheric parameter for modtarn run */
    bool verbose                /*I: value to indicate if intermediate messages 
                                     be printed */
)
{ 
    return SUCCESS;
}
