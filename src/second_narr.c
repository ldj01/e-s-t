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
int planck_eq
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
    char errstr[MAX_STR_LEN];
   
    /* Allocate memory */
    lambda = (float *)malloc(num_srs * sizeof(float));  
    if (lambda == NULL)
    {
        sprintf (errstr, "Allocating lambda memory");
        LST_ERROR (errstr, "second_narr");
    }

    /* Planck Const hecht pg, 585 ## units: Js */
    h = 6.6260755 * (float)10e-34;

    /* Boltzmann Gas Const halliday et 2001 ## units: J/K */
    k = 1.3806503 * (float)10e-23;

    /* Speed of Light ## units: m/s */
    c = 299792458.0;

    for (i = 0; i < num_srs; i++)
    {
        /* Lambda intervals of Landsat5 spectral response locations microns 
           units: m */
        lambda[i] = wavelength[i] * (float)10e-6;

        /* Compute the Planck Blackbody Eq [W/m^2 sr um] */  
        black_radiance[i] = 2.0 * h * (pow(c, 2.0)) * ((float)(10e-6) * 
            pow(lambda[i],-5.0)) * (1.0/(exp((h * c) / (lambda[i] * k * 
            temperature)) -1.0));

        /* convert to W/cm^2 sr micron to match modtran units */
        black_radiance[i] /= 10000.0; 
    }
    free(lambda);

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
void spline
(
    float *x, 
    float *y, 
    int  n, 
    float  yp1, 
    float  ypn, 
    float *y2
)
{
    int  i,k;
    float p,qn,sig,un,*u;
    char errstr[MAX_STR_LEN];

    u=(float *)malloc((unsigned) (n-1)*sizeof(float));
    if (!u)
    {
        sprintf(errstr, "Can't allocate memory");
        LST_ERROR(errstr, "second_narr");
    }

    if ((yp1 - 0.99e30) > MINSIGMA)
        y2[0]=u[0]=0.0;
    else 
    {
       	y2[0] = -0.5;
        u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
    }
    for (i=1;i<=n-2;i++) 
    {
        sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
        p=sig*y2[i-1]+2.0;
        y2[i]=(sig-1.0)/p;
        u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
        u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
    if ((ypn - 0.99e30) > MINSIGMA)
        qn=un=0.0;
    else 
    {
        qn=0.5;
        un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    }
    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
    for (k=n-2;k>=0;k--)
        y2[k]=y2[k]*y2[k+1]+u[k];
    free(u);
}

/******************************************************************************
MODULE:  splint

PURPOSE: splint uses the cubic spline generated with spline to interpolate 
         values in the XY  table

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/29/2014   Song Guo         Modified from online code
******************************************************************************/
void splint
(
    float *xa, 
    float *ya, 
    float *y2a,
    int  n, 
    float  x, 
    float *y
)
{
    int  k;
    float h,b,a;

    if ( klo < 0 )
    {
       klo=0;
       khi=n-1;
    } 
    else 
    {
        if ( x < xa[klo] ) klo=0;
        if ( x > xa[khi] ) khi=n-1;
    }
    while (khi-klo > 1) 
    {
        k=(khi+klo) >> 1;
        if ((xa[k] - x) > MINSIGMA) khi=k;
        else klo=k;
    }
    h=xa[khi]-xa[klo];
    if (h == 0.0) 
    {
        *y = 0.0;;
    } 
    else 
    {
        a=(xa[khi]-x)/h;
        b=(x-xa[klo])/h;
        *y=a*ya[klo]+b*ya[khi]+( (a*a*a-a)*y2a[klo]+
        (b*b*b-b)*y2a[khi] ) * (h*h) / 6.0;
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
    float *x,                 /*I: The tabulated X-value data */
    float *f,                 /*I: The tabulated F-value data */ 
    int nums,                 /*I: number of points */
    float *result             /*O: integraeted result */  
)
{ 
    float *temp;
    float *z;
    float xmin;
    float xmax;
    int i;
    int *ii;
    float h;
    char errstr[MAX_STR_LEN];

    /* Allocate memory */
    temp = (float *)malloc(nums * sizeof(float));  
    if (temp == NULL)
    {
        sprintf (errstr, "Allocating temp memory");
        LST_ERROR (errstr, "second_narr");
    }

    z = (float *)malloc(nums * sizeof(float));  
    if (z == NULL)
    {
        sprintf (errstr, "Allocating z memory");
        LST_ERROR (errstr, "second_narr");
    }

    ii = (int *)malloc(nums * sizeof(int));  
    if (ii == NULL)
    {
        sprintf (errstr, "Allocating ii memory");
        LST_ERROR (errstr, "second_narr");
    }

    /* integrate spectral response over wavelength */
    /* Call spline to get second derivatives, */
    spline(x, f, nums, 2.0, 2.0, temp);
    
    /* Call splint for interpolations. one-based arrays are considered */
    for (i = 0; i < nums; i++) 
    {
        splint(x, f, temp, nums, i, &z[i]);
    }

    xmin = x[0];
    xmax = x[nums-1];  
    h = (xmax - xmin) / (float) nums;

    *result = 0.0;
    /* Compute the integral using the 5-point Newton-Cotes formula */
    for (i = 0; i < nums; i++) 
    {
        ii[i] = 4 * (int)(((nums - 1) / 4) + 1);
        *result += 2.0 * h * (7.0 * (z[ii[i]-4] + z[ii[i]]) + 
                 32.0 * (z[ii[i]-3] + z[ii[i]-1]) + 12.0 * z[ii[i]-2]) / 45.0;
    }

    free(temp);
    free(z);
    free(ii);

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
    float temperature,          /*I: temperature */
    float **spectral_response,  /*I: spectral response function */
    int num_srs,                /*I: number of spectral response points */
    float *radiance             /*O: blackbody radiance */  
)
{ 
    int i;
    float rs_integral;
    float temp_integral;
    float *blackbody_radiance;
    float *product;
    char errstr[MAX_STR_LEN];
    int status;

    /* Allocate memory */
    blackbody_radiance = (float *)malloc(num_srs * sizeof(float));  
    if (blackbody_radiance == NULL)
    {
        sprintf (errstr, "Allocating blackbody_radiance memory");
        LST_ERROR (errstr, "second_narr");
    }

    product = (float *)malloc(num_srs * sizeof(float));  
    if (product == NULL)
    {
        sprintf (errstr, "Allocating product memory");
        LST_ERROR (errstr, "second_narr");
    }

    /* integrate spectral response over wavelength */
    status = int_tabulated(spectral_response[0], spectral_response[1], num_srs, 
                           &rs_integral);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Calling int_tabulated\n");
        LST_ERROR (errstr, "second_narr");              
    }

    /* using planck's equaiton to calculate radiance at each wavelength for 
       current temp*/
    status = planck_eq(spectral_response[0], num_srs, temperature, 
                       blackbody_radiance);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Calling planck_eq\n");
        LST_ERROR (errstr, "second_narr");              
    }

    /* multiply the caluclated planck radiance by the spectral reponse and integrate 
       over wavelength to get one number for current temp*/
    for (i = 0; i < num_srs; i++) 
    {
        product[i] = blackbody_radiance[i] * spectral_response[1][i];
    }

    int_tabulated(spectral_response[0], product, num_srs, &temp_integral);

    /* divide above result by integral of spectral response function */
    *radiance = temp_integral / rs_integral;

    /* Free allocated memory */
    free(blackbody_radiance);
    free(product);

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
******************************************************************************/
void interpol
(
    float **v_x,         /*I: The input vector data */     
    float *u,            /*I: The output grid points */   
    int nums,            /*I: number of input data */   
    int index,           /*I: Which index column data be used */                 
    float *r             /*O: Interpolated results */            
)
{ 
    int i;
    int ix = 0;
    float s1;
    int d;
    int m;

    m = nums - 2;

    if ((v_x[1][0] - v_x[0][0]) >= MINSIGMA)
        s1 = 1;
    else
        s1 = -1;

    for (i = 0; i < nums; i++)
    {
        d = (int) (s1 * (u[i] - v_x[ix][0]));
        if (d == 0)
            r[i] = v_x[ix][index];
        else
        {
            if (d > 0) 
            {
                while ((s1*(u[i]-v_x[ix+1][0]) > MINSIGMA) && (ix < m -2))
                    ix++;    
            }
            else
            {
                while ((s1*(u[i]-v_x[ix][0]) < MINSIGMA) && (ix > 0))
                    ix--;
            }
            r[i] = v_x[ix][index] + (u[i]-v_x[ix][0])*(v_x[ix+1][index]-v_x[ix][index]) / 
                   (v_x[ix+1][0]-v_x[ix][0]);
        }
    }

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
9/29/2014   Song Guo         Original Development
******************************************************************************/
int calculate_lobs
(
    float **modtran,            /*I: modtran results with wavelengths */
    float **spectral_response,  /*I: spectral response function */
    int num_srs,                /*I: number of spectral response points */
    int index,                  /*I: column index for data be used */
    float *radiance             /*O: LOB outputs */  
)
{ 
    int i;
    float *temp_rad;
    float rs_integral;
    float temp_integral;
    float *product;
    char errstr[MAX_STR_LEN];

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
    int_tabulated(spectral_response[0], spectral_response[1], num_srs, &rs_integral);

    /* using planck's equaiton to calculate radiance at each wavelength for 
       current temp*/
    interpol(modtran, spectral_response[0], num_srs, index, temp_rad);

    /* multiply the caluclated radiance by the spectral reponse and integrate 
       over wavelength to get one number for current temp*/
    for (i = 0; i < num_srs; i++) 
    {
        product[i] = temp_rad[i] * spectral_response[1][i];
    }

    int_tabulated(spectral_response[0], product, num_srs, &temp_integral);

    /* divide above result by integral of spectral response function */
    *radiance = temp_integral / rs_integral;

    /* Free allocated memory */
    free(temp_rad);
    free(product);

    return SUCCESS;

}

/******************************************************************************
MODULE:  get_lat_lon_height 

PURPOSE: Get latitude, longitude, and height from directory name

RETURN: None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
12/3/2014   Song Guo         Original Development
******************************************************************************/
void get_lat_lon_height 
(
    char *full_path,       /* I: Full path in case_list */
    float *lat,            /* O: latitude */
    float *lon,            /* O: longitude */
    float *height          /* O: height */
)
{
    char directory[MAX_STR_LEN];           
    char temp_value[MAX_STR_LEN];          
    char directory2[MAX_STR_LEN];            
    char full_path_copy[MAX_STR_LEN];  
    char *ptr;                  /* String pointer */

    /* Make a local copy of filename so it is not destroyed */
    strcpy(full_path_copy, full_path);

    /* Find ending '_' */
    ptr = (char *) strrchr(full_path_copy, '_');
    if (ptr != NULL)
    {
        strcpy (directory, full_path_copy);
        ptr = (char *) strrchr(directory, '_');
        ptr++;
        strcpy (full_path_copy, ptr);
        *ptr = '\0';
    }
    else
        strcpy (directory, "");

    /* Get latitude */
    ptr = (char *) strrchr(directory, '/');
    if (ptr != NULL)
    {
        ptr++;
        strncpy (temp_value, ptr, strlen(ptr)-1);
        *ptr = '\0';
    }
    *lat = atof(temp_value);

    /* Get longitude */
    ptr = (char *) strchr(full_path_copy, '/');
    if (ptr != NULL)
    {
        *(ptr++) = '\0';
        strcpy (temp_value, full_path_copy);
        strcpy (directory2, ptr);
    }
    *lon = atof(temp_value);

    /* Get height */
    ptr = (char *) strchr(directory2, '/');
    if (ptr != NULL)
    {
        *(ptr++) = '\0';
        strcpy (temp_value, directory2);
    }
    *height = atof(temp_value);

    return;
}

/******************************************************************************
MODULE:  second_narr

PURPOSE: Generate transmission, upwelled radiance, and downwelled radiance at 
         each height at each NARR point

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/29/2014   Song Guo         Original Development
******************************************************************************/
int second_narr
(
    Input_t *input,             /*I: input structure */
    int num_points,             /*I: number of narr points */
    float alb,                  /*I: albedo */ 
    char **case_list,           /*I: modtran run list */
    float **results,            /*O: atmospheric parameter for modtarn run */
    bool verbose                /*I: value to indicate if intermediate messages 
                                     be printed */
)
{ 
    FILE *fd;
    int i, j, k, m;
    float **spectral_response;
    float temp_radiance_0;
    float obs_radiance_0;
    float temp_radiance_273;
    float temp_radiance_300;
    int counter = 0;
    int place = 0;
    char command[MAX_STR_LEN];
    int index;
    float lat, lon;
    float height;
    char num_entry[MAX_STR_LEN];
    int num_entries;
    char current_file[MAX_STR_LEN];
    float **temp1;
    char zero_tape6[MAX_STR_LEN];
    char zero_temp6[MAX_STR_LEN];
    float zero_temp;
    float **current_data;
    float x[2];
    float y[2];
    float tau, lu, ld;
    float ems = 1 - alb;
    char errstr[MAX_STR_LEN];
    char *path = NULL;
    char full_path[MAX_STR_LEN];
    int status;

    path = getenv("LST_DATA");
    if (path == NULL)
    {
        sprintf (errstr, "LST_DATA environment variable is not set");
        LST_ERROR(errstr, "second_narr");
    }

    if (input->meta.inst == INST_TM && input->meta.sat == SAT_LANDSAT_5)
    {
        /* Dynamic allocate the 2d memory */
        spectral_response = (float **)allocate_2d_array(2, 171,sizeof(float)); 
        if (spectral_response == NULL)
        {
            sprintf (errstr, "Allocating spectral_response memory");
            LST_ERROR (errstr, "second_narr");
        }

        sprintf(full_path,"%s/%s", path, "L5v2.rsp");
        fd = fopen(full_path, "r");
        if (fd == NULL)
        {
            sprintf(errstr, "Can't open L5v2.rsp file");;
            LST_ERROR(errstr, "second_narr");
        }

        for (i = 0; i < 171; i++)
        {
            if (fscanf(fd, "%f %f", &spectral_response[0][i], 
                   &spectral_response[1][i]) == EOF)  
            {
                sprintf(errstr, "End of file (EOF) is met before 171 lines");;
                LST_ERROR(errstr, "second_narr");
            }
        }
        fclose(fd);
    }
    else if (input->meta.inst == INST_ETM && input->meta.sat == SAT_LANDSAT_7)
    {
        /* Dynamic allocate the 2d memory */
        spectral_response = (float **)allocate_2d_array(2, 47,sizeof(float)); 
        if (spectral_response == NULL)
        {
            sprintf (errstr, "Allocating spectral_response memory");
            LST_ERROR (errstr, "second_narr");
        }

        sprintf(full_path,"%s/%s", path, "L7.rsp");
        fd = fopen(full_path, "r");
        if (fd == NULL)
        {
            sprintf(errstr, "Can't open L7.rsp file");;
            LST_ERROR(errstr, "second_narr");
        }

        for (i = 0; i < 47; i++)
        {
            if (fscanf(fd, "%f %f", &spectral_response[0][i], 
                   &spectral_response[1][i]) == EOF)  
            {
                sprintf(errstr, "End of file (EOF) is met before 47 lines");
                LST_ERROR(errstr, "second_narr");
            }
        }
        fclose(fd);
    }
    else if (input->meta.inst == INST_OLI_TIRS && input->meta.sat == SAT_LANDSAT_8)
    {
        /* Dynamic allocate the 2d memory */
        spectral_response = (float **)allocate_2d_array(2, 101,sizeof(float)); 
        if (spectral_response == NULL)
        {
            sprintf (errstr, "Allocating spectral_response memory");
            LST_ERROR (errstr, "second_narr");
        }

        sprintf(full_path,"%s/%s", path, "L8_B10.rsp");
        fd = fopen(full_path, "r");
        if (fd == NULL)
        {
            sprintf (errstr, "Can't open L8_B10.rsp file");
            LST_ERROR(errstr, "second_narr");
        }

        for (i = 0; i < 47; i++)
        {
            if (fscanf(fd, "%f %f", &spectral_response[0][i], 
                   &spectral_response[1][i]) == EOF)  
            {
                sprintf (errstr, "End of file (EOF) is met before 101 lines");
                LST_ERROR(errstr, "second_narr");
            }
        }
        fclose(fd);
    }
    else
    {
        sprintf (errstr, "invalid instrument type");;
        LST_ERROR(errstr, "second_narr");
    }
#if 0
    /* Allocate memory */
    case_list = (char **)allocate_2d_array(MAX_STR_LEN, num_cases, sizeof(char));  
    if (case_list == NULL)
    {
        sprintf (errstr, "Allocating case_list memory");
        LST_ERROR (errstr, "second_narr");
    }

    /* Read caseList into 2d char array */
    fd = fopen("caseList", "r"); 
    if (fd == NULL)
    {
        sprintf (errstr, "Opening file: caseList\n");
        LST_ERROR (errstr, "second_narr");
    }

    /* Write out the caseList file */
    for (k = 0; k < num_poimts; k++)
    {
        fscnff(fd, "%s\n", case_list[k]);
    }

    /* Close the caseList file */
    status = fclose(fd);
    if ( status )
    {
        sprintf (errstr, "Closing file: caseList\n");
        LST_ERROR (errstr, "second_narr");
    }
#endif
    /* calculate Lt for each temperature */
    if (input->meta.inst == INST_TM && input->meta.sat == SAT_LANDSAT_5)
    {
        status = calculate_lt(273, spectral_response, 171, &temp_radiance_273);
        if (status != SUCCESS)
        { 
            sprintf (errstr, "Calling calculate_lt for 273K");
            LST_ERROR (errstr, "second_narr");
        }
        status = calculate_lt(300, spectral_response, 171, &temp_radiance_300);
        if (status != SUCCESS)
        { 
            sprintf (errstr, "Calling calculate_lt for 300K");
            LST_ERROR (errstr, "second_narr");
        }
    }
    else if (input->meta.inst == INST_ETM && input->meta.sat == SAT_LANDSAT_7)
    {
        status = calculate_lt(273, spectral_response, 47, &temp_radiance_273);
        if (status != SUCCESS)
        { 
            sprintf (errstr, "Calling calculate_lt for 273K");
            LST_ERROR (errstr, "second_narr");
        }
        status = calculate_lt(300, spectral_response, 47, &temp_radiance_300);
        if (status != SUCCESS)
        { 
            sprintf (errstr, "Calling calculate_lt for 300K");
            LST_ERROR (errstr, "second_narr");
        }
    }
    else if (input->meta.inst == INST_OLI_TIRS && input->meta.sat == SAT_LANDSAT_8)
    {
        status = calculate_lt(273, spectral_response, 101, &temp_radiance_273);
        if (status != SUCCESS)
        { 
            sprintf (errstr, "Calling calculate_lt for 273K");
            LST_ERROR (errstr, "second_narr");
        }
        status = calculate_lt(300, spectral_response, 101, &temp_radiance_300);
        if (status != SUCCESS)
        { 
            sprintf (errstr, "Calling calculate_lt for 300K");
            LST_ERROR (errstr, "second_narr");
        }
    }
    else
    {
        sprintf (errstr, "invalid instrument type");
        LST_ERROR(errstr, "second_narr");
    }

    /* iterate through all points in the scene and all heights in one point */
    for (i = 0; i < num_points; i++)
    {
        for (j = 0; j < NUM_ELEVATIONS; j++)
        {
            /* determine current latlon and height depends on number of steps 
               in path */
            get_lat_lon_height(case_list[counter], &lat, &lon, &height);

#if 0
            sprintf(command, "echo %s|tr '/' '\\n'",case_list[counter]);
            printf("command=%s\n",command);

            /* Open the command for reading. */
            fd = popen(command, "r");
            if (fd == NULL) 
            {
                sprintf (errstr, "Allocating results memory");
                LST_ERROR (errstr, "second_narr");
            }

            index = 0;
            /* Read the output a line at a time - output it. */
            while (fgets(spec[index], sizeof(spec[index])-1, fd) != NULL) 
                index++;

            /* close the command file */
            pclose(fd);

            /* Find the right location for lat_lon in the spec */

            /* Split lat and lon */
            sprintf(command, "echo %s|tr '_' '\\n'",spec[4]);

            /* Open the command for reading. */
            fd = popen(command, "r");
            if (fd == NULL) 
            {
                sprintf (errstr, "Allocating results memory");
                LST_ERROR (errstr, "second_narr");
            }

            index = 0;
            /* Read the output a line at a time - output it. */
            while (fgets(coordinates[index], sizeof(coordinates[index])-1, fd) 
                   != NULL) 
                index++;

            /* close the command file */
            pclose(fd);

            /* get lat, lon and height */
            sscanf(coordinates[0], "%f", &lat);          
            sscanf(coordinates[1], "%f", &lon);          
            sscanf(spec[5], "%f", &height);          
#endif
            /* determine number of entries in current file */
            sprintf(command, "wc %s/parsed | awk '{print $1}'", case_list[counter]);

            /* Open the command for reading. */
            fd = popen(command, "r");
            if (fd == NULL) 
            {
                sprintf (errstr, "Getting number of entries");
                LST_ERROR (errstr, "second_narr");
            }

            /* Read the output a line at a time - output it. */
            fgets(num_entry, sizeof(num_entry)-1, fd); 

            /* close the command file */
            pclose(fd);

            /* get number of entries */
            sscanf(num_entry, "%d", &num_entries);  

            printf("i,j,num_entry,num_entries=%d,%d,%s,%d\n",
                    i,j,num_entry,num_entries);

            /* for each height, read in radiance inforomation for three modtran 
               runs, columns of array are organized:
               wavelength | 273,0.0 | 310,0.0 | 000,0.1 */
            current_data = (float **)allocate_2d_array(num_entries, 4, sizeof(float)); 
            if (current_data == NULL)
            {
                sprintf (errstr, "Allocating current_data memory");
                LST_ERROR (errstr, "second_narr");
            }

            temp1 = (float **)allocate_2d_array(num_entries, 2, sizeof(float)); 
            if (temp1 == NULL)
            {
                sprintf (errstr, "Allocating temp1 memory");
                LST_ERROR (errstr, "second_narr");
            }

            index = 0;
            /* iterate through three pairs of parameters */
            for (k = 0; k < 2; k++)
            {
                /* define current file */
                sprintf(current_file, "%s/parsed", case_list[counter]);

                fd = fopen(current_file, "r");
                if (fd == NULL)
                {
                    sprintf (errstr, "Can't open current_file file");
                    LST_ERROR(errstr, "second_narr");
                }

                for (m = 0; m < num_entries ; m++)
                {
                    if (fscanf(fd, "%f %f", &temp1[m][0], &temp1[m][1]) == EOF)  
                    {
                        sprintf (errstr, "End of file (EOF) is met before num_entries "
                                 " lines");
                        LST_ERROR(errstr, "second_narr");
                    }
                }
                fclose(fd);

                /* put arrays into data array for current point at current height */
                if (index == 0) 
                {
                    for (m = 0; m < num_entries ; m++)
                    {
                        current_data[m][0] = temp1[m][0];
                    }
                    index++;
                }

                for (m = 0; m < num_entries ; m++)
                {
                    current_data[m][index] = temp1[m][1];
                }

                if (k == 2)
                {
                    /* determine temperature at lowest atmospheric layer 
                       (when MODTRAN is run at 0K) */
                    sprintf(zero_tape6, "%s/tape6", case_list[counter]);
                    sprintf(command, "grep \"TARGET-PIXEL (H2) SURFACE TEMPERATURE\" %s | " 
                                     "awk '{print $7}'", zero_tape6 );

                    /* Open the command for reading. */
                    fd = popen(command, "r");
                    if (fd == NULL) 
                    {
                        sprintf (errstr, "Popening command");
                        LST_ERROR (errstr, "second_narr");
                    }

                    /* Read the output a line at a time - output it. */
                    while (fgets(zero_temp6, sizeof(zero_temp6)-1, fd) != NULL) 

                    /* close the command file */
                    pclose(fd);

                    /* get zero_temp */
                    sscanf(zero_temp6, "%f", &zero_temp);          
                }
                counter++;
                index++;
            }

            /* parameters from 3 modtran runs
               Lobs = Lt*tau + Lu; m = tau; b = Lu; */
            x[0] = temp_radiance_273;
            x[1] = temp_radiance_300;
            status = calculate_lobs(current_data, spectral_response, 
                                    num_entries, 1, &y[0]);
            if (status != SUCCESS)
            {  
                sprintf (errstr, "Calling calculate_lob 1");
                LST_ERROR (errstr, "second_narr");
            }
            status = calculate_lobs(current_data, spectral_response, 
                                    num_entries, 2, &y[1]);
            if (status != SUCCESS)
            {  
                sprintf (errstr, "Calling calculate_lob 2");
                LST_ERROR (errstr, "second_narr");
            }         

            /* Free the allocated memory in the loop*/
            status = free_2d_array((void **)current_data);
            if (status != SUCCESS)
            {
                sprintf (errstr, "Freeing memory: current_data\n");
                LST_ERROR (errstr, "second_narr");              
            }

            status = free_2d_array((void **)temp1);
            if (status != SUCCESS)
            {
                sprintf (errstr, "Freeing memory: temp\n");
                LST_ERROR (errstr, "second_narr");              
            }

            /* Implement a = INVERT(TRANSPOSE(x)##x)##TRANSPOSE(x)##y 
               Note: I slove the two equations analytically */
            tau = (y[0] - y[1]) / (x[0] - x[1]);
            lu = (x[1]*y[0] - x[0]*y[1]) / (x[0] - x[1]);

            /* determine Lobs and Lt when modtran was run a 0 K - calculate downwelled */
            if (input->meta.inst == INST_TM && input->meta.sat == SAT_LANDSAT_5)
            {
                status = calculate_lt(zero_temp, spectral_response, 171, &temp_radiance_0);
                if (status != SUCCESS)
                { 
                    sprintf (errstr, "Calling calculate_lt for 0K");
                    LST_ERROR (errstr, "second_narr");
                }
            }
            else if (input->meta.inst == INST_ETM && input->meta.sat == SAT_LANDSAT_7)
            {
                status = calculate_lt(zero_temp, spectral_response, 47, &temp_radiance_0);
                if (status != SUCCESS)
                { 
                    sprintf (errstr, "Calling calculate_lt for 0K");
                    LST_ERROR (errstr, "second_narr");
                }
            }
            else if (input->meta.inst == INST_OLI_TIRS && input->meta.sat == SAT_LANDSAT_8)
            {
                status = calculate_lt(zero_temp, spectral_response, 101, &temp_radiance_0);
                if (status != SUCCESS)
                { 
                    sprintf (errstr, "Calling calculate_lt for 0K");
                    LST_ERROR (errstr, "second_narr");
                }
            }
            else
            {
                sprintf (errstr, "invalid instrument type");
                LST_ERROR(errstr, "second_narr");
            }
            status = calculate_lobs(current_data, spectral_response, 
                                    num_entries, 3, &obs_radiance_0);
            if (status != SUCCESS)
            {  
                sprintf (errstr, "Calling calculate_lob 2");
                LST_ERROR (errstr, "second_narr");
            }         

            /* Ld = (((Lobs-Lu)/tau) - (Lt*ems))/(1-ems) */
            ld = (((obs_radiance_0-lu) / tau) - (temp_radiance_0 * ems)) / (1 - ems); 

            /* put results into results array */
            results[i * NUM_ELEVATIONS + j][0] = lat;
            results[i * NUM_ELEVATIONS + j][1] = lon; 
            results[i * NUM_ELEVATIONS + j][2] = height;
            results[i * NUM_ELEVATIONS + j][3] = tau;
            results[i * NUM_ELEVATIONS + j][4] = lu;
            results[i * NUM_ELEVATIONS + j][5] = ld;
            place++;       
        }
    }
#if 0
    /* Free allocated memory */
    status = free_2d_array((void **)spectral_response);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: spectral_response\n");
        LST_ERROR (errstr, "second_narr");              
    }

    /* write results to a file */
    sprintf(full_path,"%s/%s", path, "atmosphericParameters.txt");
    fd = fopen(full_path, "w");
    if (fd == NULL)
    {
        sprintf (errstr, "Can't open atmosphericParameters.txt file");
        LST_ERROR(error_string, "second_narr");
    }
    for (k = 0; k < num_points * NUM_ELEVATIONS; k++)
    {
        fprintf(fd, "%f,%f,%f,%f,%f,%f\n", 
                results[0][k], results[1][k], results[2][k], 
                results[3][k], results[4][k], results[5][k]);
    }
    fclose(fd);
#endif
    return SUCCESS;
}
