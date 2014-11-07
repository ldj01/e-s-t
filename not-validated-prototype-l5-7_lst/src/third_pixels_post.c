#include <string.h>
#include <stdarg.h>
#include <math.h>
#include "2d_array.h"
#include "input.h"
#include "const.h"
#include "scene_based_lst.h"

#ifndef max
#define max(a,b) (((a) (b)) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

int partition(float arr[], int left, int right)
{
      int i = left, j = right;
      int tmp;
      int pivot = arr[(left + right) / 2];

      while (i <= j) {
            while (arr[i] < pivot)
                  i++;
            while (arr[j] > pivot)
                  j--;
            if (i <= j) {
                  tmp = arr[i];
                  arr[i] = arr[j];
                  arr[j] = tmp;
                  i++;
                  j--;
            }
      }

      return i;
}

void quick_sort(float arr[], int left, int right) 
{
      int index = partition(arr, left, right);

      if (left < index - 1)
       quick_sort(arr, left, index - 1);
      if (index < right)
       quick_sort(arr, index, right);
}

/******************************************************************************
MODULE:  distance_in_utm

PURPOSE: Calculate distances between UTM coordiantes

RETURN: NONE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
10/10/2014   Song Guo         Original Development
******************************************************************************/
void distance_in_utm
(
    float e1,
    float n1,
    float e2,
    float n2,
    float *dis
)
{
    float s = 0.9996;   /* scale factor */
    float r = 6378137;  /* earth radius */
    float sr1, sr2, sr3;
    float edist;

    sr1 = s / (cos(e1/r));
    sr2 = s / (cos(((e2 - e1) / 6.0) / r));
    sr3 = s / (cos(e2 / r));
   
    edist = ((e2-e1)/6.0)*(sr1+4*sr2+sr3);
   
    *dis = sqrt(edist*edist+(n2-n1)*(n2-n1));
}

/******************************************************************************
MODULE:  convert_ll_utm

PURPOSE: convert digital counts to radiance for thermal band
         [unit: W m^(-2) sr^(-1) mu^(-1) ]

RETURN: NONE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/30/2014   Song Guo         Original Development
******************************************************************************/
void convert_ll_utm
(
    Input_t *input,
    float *lat,
    float *lon,
    int num_points,
    float **narr_utm
)
{
    int i;
    float a = 6378137.0;    /* equatorial radius */
    float b = 6356752.3;    /* polar radius */
    float k0 = 0.9996;      /* scale factor */
    float e;                /* eccentricity */
    float eprimesqrd, n, rho, nu;
    float a0, b0, c0, d0, e0;
    float ki, kii, kiii, kiv, kv;
    float zone_cm, zone_cm_rad;
    float delta_lon;
    float p, lat_rad, lon_rad;
    float s;
    float northing, easting; 

    /* calculate zone central meridian in degrees and radians */
    zone_cm = (float) (6 * input->meta.zone - 183);
    zone_cm_rad =  zone_cm * PI / 180.0;

    e = sqrt(1.0-pow(b/a, 2));
    eprimesqrd = e*e/(1.0-e*e);
    n = (a-b)/(a+b);

    /* calculate meridional arc length */   
    a0 = a*(1.0-n+(5.0*n*n/4.0)*(1.0-n)+(81.0* (float)pow(n,4)/64.0)*(1.0-n));
    b0 = (3.0*a*n/2.0)*(1.0-n-(7.0*n*n/8.0)*(1.0-n)+55.0*(float)pow(n,4)/64.0);
    c0 = (15.0*a*n*n/16.0)*(1.0-n+(3.0*n*n/4.0)*(1.0-n));
    d0 = (35*a*(float)pow(n,3)/48.0)*(1.0-n+11.0*n*n/16.0);
    e0 = (315.0*a*(float)pow(n,4)/51.0)*(1.0-n);

    for ( i = 0; i < num_points; i++)
    {   
        delta_lon = lon[i] - zone_cm;
        p = delta_lon * PI/ 180.0;
        /* convert lat and lon points from decimal degrees to radians */
        lat_rad = lat[i]*PI/180.0;
        lon_rad = lon[i]*PI/180.0;

        rho = a*(1.0-e*e)/(pow((1.0-(e*sin(lat_rad))*(e*sin(lat_rad))),
              (3.0/2.0)));
        nu = a/pow((1.0-(e*sin(lat_rad))*(e*sin(lat_rad))),(1.0/2.0));


        s = a0*lat_rad - b0*sin(2*lat_rad) + c0*sin(4*lat_rad) - 
            d0*sin(6*lat_rad) + e0*sin(8*lat_rad);

        /* coefficients for UTM coordinates */   
        ki = s*k0;
        kii = nu*sin(lat_rad)*cos(lat_rad)*k0/2.0;
        kiii = (pow(nu*sin(lat_rad)*cos(lat_rad),3)/24.0)*
            pow((5-tan(lat_rad)),2)+9.0*eprimesqrd*pow(cos(lat_rad),2)+
            4.0*eprimesqrd*eprimesqrd*pow(cos(lat_rad),4)*k0;
        kiv = nu*cos(lat_rad)*k0;
        kv = pow(cos(lat_rad),3)*(nu/6.0)*(1-tan(lat_rad)*tan(lat_rad)+
             eprimesqrd*cos(lat_rad)*cos(lat_rad))*k0;   

        /* calculate UTM coordinates */   
        northing = (ki+kii*p*p+kiii*pow(p,4));
        easting = 500000.0 + (kiv*p+kv*pow(p,3));

        /* assign to narr_utm array*/ 
        narr_utm[0][i] = easting;
        narr_utm[0][i] = northing;
    }
}

/******************************************************************************
MODULE:  dn_to_radiance

PURPOSE: convert digital counts to radiance for thermal band
         [unit: W m^(-2) sr^(-1) mu^(-1) ]

RETURN: NONE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/30/2014   Song Guo         Original Development
******************************************************************************/
void dn_to_radiance
(
    Input_t *input
)
{
    int pix;
    for (pix = 0; pix < input->size_th.s; pix++)
    { 
        if (input->therm_buf[pix] == 0)
            input->therm_buf[pix] = 0;
        else
        {
            input->therm_buf[pix] = (int16) (input->meta.gain_th * 
                   (float)input->therm_buf[pix] + 
                   input->meta.bias_th);
        }   
        if (input->meta.inst == INST_TM && input->meta.sat == SAT_LANDSAT_5)
         input->therm_buf[pix] = (int16) (input->therm_buf[pix] + 0.044);
    }
}

/******************************************************************************
MODULE:  interpolate_to_height

PURPOSE: Inteprolate to height of current pixel

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
10/29/2014   Song Guo         
******************************************************************************/
void interpolate_to_height
(
    int current_index,
    float *height,              /*I: list of height of one location */
    float *atm1,                /*I: atmospheric parameter colum1 */
    float *atm2,                /*I: atmospheric parameter colum2 */
    float *atm3,                /*I: atmospheric parameter colum3 */
    float interpolate_to,       /*I: current landsat pixel height */
    float *at_height            /*O: interpolated height */  
)
{ 
    int i;
    int below = 0;
    int above = 0;
    int close_below;
    int close_above;
    float under_height;
    float under_variables[3];
    float over_height;
    float over_variables[3];
    float m[3];
    float b[3];

    /* Determine points below/above interpolate_to */
    for (i = current_index; i < current_index + NUM_ELEVATIONS - 1; i++)
    {
        if ((height[i] - interpolate_to) < MINSIGMA)
            below++;
        else
            above++;
    }

    if (below == 0)
        close_below = 0;
    else
        close_below = below;
    under_height = height[close_below];
    under_variables[0] = atm1[close_below];
    under_variables[1] = atm2[close_below];
    under_variables[2] = atm3[close_below];

    if (above == 0)
        close_above = NUM_ELEVATIONS - 1;
    else
        close_below = above;
    over_height = height[close_below];
    over_variables[0] = atm1[close_below];
    over_variables[1] = atm2[close_below];
    over_variables[2] = atm3[close_below];

    for (i = 0; i < 3; i++)
    {
        if (close_above == close_below)
            at_height[i] = under_variables[i];
        else
        {
            m[i] = (over_variables[i] - under_variables[i]) / 
                (over_height - under_height);
            b[i] = over_variables[i] - m[i] * over_height;
            at_height[i] = m[i] * interpolate_to + b[i];
        }
    }
}

/******************************************************************************
MODULE:  interpolate_to_location

PURPOSE: Inteprolate to location of current pixel

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
10/29/2014   Song Guo         
******************************************************************************/
void interpolate_to_location
(
    float **coordinates,        /*I: current pixel coordinates */
    float **at_height,          /*I: current height atmospheric results */
    float interpolate_easting,  /*I: interpolate to easting */
    float interpolate_northing, /*I: interpolate to northing */
    float *parameters           /*O: interpolated pixel atmospheric parameters */
)
{ 
    int i, j;
    float h[4];
    float w[4];
    float total = 0.0;

    /* shepard's method */
    for (i = 0; i < 4; i++)
    {
        h[i] = (coordinates[i][0] - interpolate_easting) * (coordinates[i][0] - 
            interpolate_easting) + (coordinates[i][1] - interpolate_northing) * 
            (coordinates[i][1] - interpolate_northing);
    }
    quick_sort(h, 0, 3);

    for (i = 0; i < 4; i++)
    {
        total += 1.0 / h[i];
    }

    for (i = 0; i < 4; i++)
    {
        w[i] = (1.0 / h[i]) / total;
    }

    for (i = 0; i < 3; i++)
    {
        parameters[i] = 0.0;
        for (j = 0; j < 4; j++) 
        {
            parameters[i] += w[j] * at_height[i][j];
        }
    } 
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
10/8/2014   Song Guo         Original Development
******************************************************************************/
int third_pixels_post
(
    Input_t *input,             /*I: input structure */
    int num_points,             /*I: number of narr points */
    char *dem_infile,           /*I: address of input DEM filename */
    char *emi_infile,           /*I: address of input Emissivity filename */
    float **modtran_results,    /*I: atmospheric parameter for modtarn run */
    bool verbose                /*I: value to indicate if intermediate messages 
                                     be printed */
)
{ 
    int row, col;
    int16 *dem=NULL;         /* input DEM data (meters) */
    int offset;              /* offset in the raw binary DEM file to seek to
                                to begin reading the window in the DEM */
    FILE *dem_fptr=NULL;     /* input scene-based DEM file pointer */
    int **eye;
    int **jay;
    float **lat;
    float **lon;
    float narr_ul_lat;
    float narr_ul_lon;
    float narr_lr_lat;
    float narr_lr_lon;
    int in_counter = 0;
    int counter;
    int max_eye;
    int min_eye;
    int max_jay;
    int min_jay;
    int num_eyes;
    int num_jays;
    float *narr_lat;
    float *narr_lon;
    float **narr_utm;
    float **east_grid;
    float **north_grid;
    int first_line;
    float current_easting;
    float current_northing;
    float *distances; 
    float *distance_copy; 
    int i, j, g, k, n;
    int closest[6];
    float easting_near[6];
    float northing_near[6];
    int below[6];
    int min_inlat, min_inlon;
    float **coordinates;
    int indices[4];
    float stay_up, stay_down, stay_right; 
    float move_up, move_down, move_right;
    float **at_height; 
    int current_index;
    float current_location[6][NUM_ELEVATIONS];
    float parameters[3];
    char therm_fname[] = "therm_radiance";
    char up_fname[] = "upwell_radiance";
    char down_fname[] = "downwell_radiance";
    char trans_fname[] = "atmos_transmittance";
    FILE *therm_fptr = NULL;
    FILE *trans_fptr = NULL;
    FILE *up_fptr = NULL;
    FILE *down_fptr = NULL;
    FILE *fd = NULL;
    float **landsat_results;
    char errstr[MAX_STR_LEN];           /* error string */
    char full_path[MAX_STR_LEN];         
    char *lst_path = NULL;     
    int status;    

    /* Dynamic allocate the 2d memory */
    eye = (int **)allocate_2d_array(NARR_ROW, NARR_COL, sizeof(int)); 
    if (eye == NULL)
    {
        sprintf (errstr, "Allocating eye memory");
        LST_ERROR (errstr, "third_pixels_post");
    }

    jay = (int **)allocate_2d_array(NARR_ROW, NARR_COL, sizeof(int)); 
    if (jay == NULL)
    {
        sprintf (errstr, "Allocating jay memory");
        LST_ERROR (errstr, "third_pixels_post");
    }

    lat = (float **)allocate_2d_array(NARR_ROW, NARR_COL, sizeof(float)); 
    if (lat == NULL)
    {
        sprintf (errstr, "Allocating lat memory");
        LST_ERROR (errstr, "third_pixels_post");
    }

    lon = (float **)allocate_2d_array(NARR_ROW, NARR_COL, sizeof(float)); 
    if (lon == NULL)
    {
        sprintf (errstr, "Allocating lon memory");
        LST_ERROR (errstr, "third_pixels_post");
    }

    lst_path = getenv("LST_DATA");
    if (lst_path == NULL)
    {
        sprintf (errstr, "LST_DATA environment variable is not set");
        LST_ERROR(errstr, "third_pixels_post");
    }

    sprintf(full_path,"%s/%s",lst_path,"coordinates.txt");
    fd = fopen(full_path, "r");
    if (fd == NULL)
    {
        sprintf (errstr, "Can't open coordinates.txt file");
        LST_ERROR(errstr, "third_pixels_post");
    }

    for (j = 0; j < NARR_COL; j++)
    {
        for (i = 0; i < NARR_ROW; i++)
        {
            if (fscanf(fd, "%d %d %f %f", &eye[i][j], &jay[i][j], &lat[i][j], 
                       &lon[i][j]) == EOF)
            {
                sprintf (errstr, "End of file (EOF) is met before NARR_ROW * "
                      "NARR_COL lines");
                LST_ERROR(errstr, "firstfile");
            }

            /* adjust longitude range to [-180, 180] */ 
            if ((lon[i][j] - 180.0) > MINSIGMA)
                lon[i][j] = 360.0 - lon[i][j];
            else
                lon[i][j] = - lon[i][j];
        }
    }
    fclose(fd);

    /* expand range to include NARR points outside image for edge pixels */
    narr_ul_lat = input->meta.ul_geo_corner.lat + 0.5;
    narr_ul_lon = input->meta.ul_geo_corner.lon - 0.5;
    narr_lr_lat = input->meta.lr_geo_corner.lat - 0.5;
    narr_lr_lon = input->meta.lr_geo_corner.lon + 0.5;

    /* determine what points in the NARR dataset fall within the Landsat image using 
       logical operators lessThanLat and greaterThanLat are values where the NARR values 
       are less than or greater than the edges of the Landsat corners values respectively
       pixels that are true in both fall within the Landsat scene
       the same thing is done with longitude values */
    max_eye = 0;
    min_eye = 1000;
    max_jay = 0;
    min_jay = 1000;
    for (i = 0; i < NARR_ROW - 1; i++)
    {
        for (j = 0; j < NARR_COL - 1; j++)
        {
            if ((lat[i][j] - narr_ul_lat) < MINSIGMA && 
                (lat[i][j] - narr_lr_lat) > MINSIGMA &&
                (lon[i][j] - narr_lr_lon) < MINSIGMA && 
                (lon[i][j] - narr_ul_lon) > MINSIGMA)
            {
                max_eye = max(max_eye, eye[i][j]); 
                min_eye = min(min_eye, eye[i][j]); 
                max_jay = max(max_jay, jay[i][j]); 
                min_jay = min(min_jay, jay[i][j]); 
                in_counter++;
            }
        }
    }    
    max_eye--;
    min_eye--;
    max_jay--;
    min_jay--;
    num_eyes = max_eye - min_eye + 1;
    num_jays = max_jay - min_jay + 1;

    /* Allocate memory for height of NARR points within the rectangular */
    narr_lat = (float *)malloc(num_points * sizeof(float)); 
    if (narr_lat == NULL)
    {
        sprintf (errstr, "Allocating narr_lat memory");
        LST_ERROR (errstr, "third_pixels_post");
    }

    narr_lon = (float *)malloc(num_points * sizeof(float)); 
    if (narr_lon == NULL)
    {
        sprintf (errstr, "Allocating narr_lon memory");
        LST_ERROR (errstr, "third_pixels_post");
    }

    narr_utm = (float **)allocate_2d_array(num_points, num_points, 
               sizeof(float)); 
    if (narr_utm == NULL)
    {
        sprintf (errstr, "Allocating narr_utm memory");
        LST_ERROR (errstr, "third_pixels_post");
    }

    /* extract coordinates within the NARR rectangle */
    for (i = min_eye; i <= max_eye; i++)
    {
        for (j = min_jay; j <= max_jay; j++)
        {
            narr_lat[(i-min_eye) * num_jays + (j-min_jay)] = lat[i][j]; 
            narr_lon[(i-min_eye) * num_jays + (j-min_jay)] = lat[i][j]; 
        }
    }

    /* Release memory */
    status = free_2d_array((void **)eye);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: eye\n");
        RETURN_ERROR (errstr, "third_pixels_post", FAILURE);              
    }

    status = free_2d_array((void **)jay);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: jay\n");
        RETURN_ERROR (errstr, "third_pixels_post", FAILURE);              
    }

    status = free_2d_array((void **)lat);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: lat\n");
        RETURN_ERROR (errstr, "third_pixels_post", FAILURE);              
    }

    status = free_2d_array((void **)lon);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: lon\n");
        RETURN_ERROR (errstr, "third_pixels_post", FAILURE);              
    }

    /* Convert NARR lat/lon to UTM */
    convert_ll_utm(input, narr_lat, narr_lon, num_points, narr_utm);

    /* Dynamic allocate the 2d memory */
    east_grid = (float **)allocate_2d_array(NARR_ROW, NARR_COL, sizeof(float)); 
    if (eye == NULL)
    {
        sprintf (errstr, "Allocating easy_grid memory");
        LST_ERROR (errstr, "third_pixels_post");
    }

    north_grid = (float **)allocate_2d_array(NARR_ROW, NARR_COL, sizeof(float)); 
    if (jay == NULL)
    {
        sprintf (errstr, "Allocating north_grid memory");
        LST_ERROR (errstr, "third_pixels_post");
    }

    for (i = 0; i < NARR_ROW - 1; i++)
    {
        for (j = 0; j < NARR_COL - 1; j++)
        {
            east_grid[i][j] = narr_utm[0][i*NARR_COL+j];
            north_grid[i][j] = narr_utm[1][i*NARR_COL+j];
        }
    }

    status = free_2d_array((void **)narr_utm);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: narr_utm\n");
        RETURN_ERROR (errstr, "third_pixels_post", FAILURE);              
    }

    /* Open the DEM for reading raw binary */
    dem_fptr = fopen (dem_infile, "rb");
    if (dem_fptr == NULL)
    {
        sprintf (errstr, "Error opening the DEM file: %s", dem_infile);
	LST_ERROR (errstr, "third_pixels_post");
    }

    /* Allocate memory for the DEM */
    dem = (int16 *) calloc (input->size_th.s, sizeof(int16));
    if (dem == NULL)
    {
        sprintf (errstr, "Error allocating memory for the DEM data");
        LST_ERROR (errstr, "third_pixels_post");
    }

    /* Open the intermediate binary files for writing
       Note: Needs to be deleted before release */
    therm_fptr = fopen(therm_fname, "wb"); 
    if (therm_fptr == NULL)
    {
        sprintf(errstr, "Opening report file: %s", therm_fname);
        LST_ERROR (errstr, "third_pixels_post");
    }

    trans_fptr = fopen(trans_fname, "wb"); 
    if (trans_fptr == NULL)
    {
        sprintf(errstr, "Opening report file: %s", trans_fname);
        LST_ERROR (errstr, "third_pixels_post");
    }

    up_fptr = fopen(up_fname, "wb"); 
    if (up_fptr == NULL)
    {
        sprintf(errstr, "Opening report file: %s", up_fname);
        LST_ERROR (errstr, "third_pixels_post");
    }

    down_fptr = fopen(down_fname, "wb"); 
    if (down_fptr == NULL)
    {
        sprintf(errstr, "Opening report file: %s", down_fname);
        LST_ERROR (errstr, "third_pixels_post");
    }

    /* Allocate memory for landsat_results */
    landsat_results = (float **)allocate_2d_array(3, input->size_th.s, sizeof(float)); 
    if (landsat_results == NULL)
    {
        sprintf (errstr, "Allocating landsat_results memory");
        LST_ERROR (errstr, "third_pixels_post");
    }

    /* Allocate memory for coordinates */
    coordinates = (float **)allocate_2d_array(4, 2, sizeof(float)); 
    if (coordinates == NULL)
    {
        sprintf (errstr, "Allocating coordinates memory");
        LST_ERROR (errstr, "third_pixels_post");
    }

    /* Allocate memory for at_height */
    at_height = (float **)allocate_2d_array(4, 3, sizeof(float)); 
    if (at_height == NULL)
    {
        sprintf (errstr, "Allocating at_height memory");
        LST_ERROR (errstr, "third_pixels_post");
    }

    if (verbose)
        printf("Iterate through all rows in landsat scene\n");
    /* Loop through each line in the image */
    for (row = 0; row < input->size_th.l; row++)
    {
        /* Print status on every 1000 lines */
        if (!(row%1000)) 
        {
            if (verbose)
            {
                printf ("Processing line %d\r",row);
                fflush (stdout);
            }
        }

	/* Read the input thermal band -- data is read into input->therm_buf */
	if (!GetInputThermLine(input, row))
        {
	    sprintf (errstr, "Reading input thermal data for line %d", row);
	    LST_ERROR (errstr, "third_pixels_post");
	}
        dn_to_radiance(input);

        /* Start reading DEM from the start_line */
        offset = sizeof (int16) * row * input->size_th.s;
        fseek (dem_fptr, offset, SEEK_SET);
        if (fread (dem, sizeof (int16), input->size_th.s, dem_fptr)
            != input->size_th.s)
        {
            sprintf (errstr, "Error reading values from the DEM file "
                "starting at line %d.", row);
	    LST_ERROR (errstr, "third_pixels_post");
        } 

        /* Can also read in one line of DEM data here*/

        /* set first_line be 1 */
        first_line = 1;
        for (col = 0; col < input->size_th.s; col++)
        {
            if (input->therm_buf != 0)
            {
                /* determine UTM coordinates of current pixel */
                current_easting = input->meta.ul_map_corner.x + row * 
                    input->meta.pixel_size[0];
                current_northing = input->meta.ul_map_corner.y + row * 
                    input->meta.pixel_size[1];

                if (first_line == 1)
                {
                    /* compute distance between current pixel and each narr points 
                       in UTM coordinates 
                       Note: consider only calculating points within a small range 
                             nearby */
                    distances = malloc(num_points * sizeof(float));
                    if (distances == NULL)
                    {
                        sprintf (errstr, "Allocating jay memory");
                        LST_ERROR (errstr, "third_pixels_post");
                    }
                    distance_copy = malloc(num_points * sizeof(float));
                    if (distances == NULL)
                    {
                        sprintf (errstr, "Allocating jay memory");
                        LST_ERROR (errstr, "third_pixels_post");
                    }
                    for (g = 0; g < num_points; g++)
                    {
                        distance_in_utm(narr_utm[0][g], narr_utm[1][g], current_easting, 
                            current_northing, &distances[g]);
                        distance_copy[g] = distances[g];
                    }

                    /* find indices of 6 closet points */
                    quick_sort(distances, 0, num_points-1);

                    n = 0;
                    /* find indices of 6 closest points */
                    for (g = 0; g < 6; g++)
                    {
                        for ( k = 0; k < num_points -1; k++)
                        {
                            if (distances[g] == distance_copy[k])
                            {
                                closest[g] = k;
                                easting_near[g] = narr_utm[0][k];
                                northing_near[g] = narr_utm[1][k];
                                if ((northing_near[g] - current_northing) <= 
                                    MINSIGMA)
                                {
                                    below[n] = k; 
                                    n++;
                                }
                            }
                        }
                    }

                    min_inlat = narr_lat[closest[below[0]]];
                    min_inlon = narr_lon[closest[below[0]]];
                    for (g = 1; g < n; g++)
                    {
                        min_inlat = min(min_inlat, narr_lat[closest[below[g]]]);
                        min_inlon = min(min_inlon, narr_lon[closest[below[g]]]);
                    }

                    /* extract UTM coordinates of four points to be interpolated 
                       and build array */
                    coordinates[0][0] =  east_grid[min_inlat][min_inlon];
                    coordinates[0][1] =  north_grid[min_inlat][min_inlon];
                    coordinates[1][0] =  east_grid[min_inlat][min_inlon+1];
                    coordinates[1][1] =  north_grid[min_inlat][min_inlon+1];
                    coordinates[2][0] =  east_grid[min_inlat+1][min_inlon];
                    coordinates[2][1] =  north_grid[min_inlat+1][min_inlon];
                    coordinates[3][0] =  east_grid[min_inlat+1][min_inlon+1];
                    coordinates[3][1] =  north_grid[min_inlat+1][min_inlon+1];

                    /* determine index of four points in order to pull from 
                       atmospheric parameter file */
                    indices[0] = min_inlat * num_jays + min_inlon;
                    indices[1] = (min_inlat + 1) * num_jays + min_inlon;
                    indices[2] = min_inlat * num_jays + min_inlon + 1;
                    indices[3] = (min_inlat + 1) * num_jays + min_inlon + 1;

                    /* set firstInLine variable to false */
                    first_line = 0;
                }
                else
                {
                    /* given indices of previous pixel, there are six possible 
                       quads to move into check 6 distances to determine new 
                       upperleft corner*/
                 distance_in_utm(east_grid[min_inlat][min_inlon], 
                         north_grid[min_inlat][min_inlon], 
                         current_easting, current_northing, &stay_right);

                 if ((min_inlon + 2) < num_jays)
                     distance_in_utm(east_grid[min_inlat][min_inlon+2], 
                         north_grid[min_inlat][min_inlon+2], 
                         current_easting, current_northing, &move_right);
                 else
                     stay_right = (float) SHRT_MAX;

                 distance_in_utm(east_grid[min_inlat+1][min_inlon+1], 
                         north_grid[min_inlat+1][min_inlon+1], 
                         current_easting, current_northing, &stay_up);

                 distance_in_utm(east_grid[min_inlat-1][min_inlon+1], 
                         north_grid[min_inlat-1][min_inlon+1], 
                         current_easting, current_northing, &move_up);

                 distance_in_utm(east_grid[min_inlat][min_inlon+1], 
                         north_grid[min_inlat][min_inlon+1], 
                         current_easting, current_northing, &stay_down);

                 if ((min_inlat + 2) < num_eyes)
                     distance_in_utm(east_grid[min_inlat+2][min_inlon], 
                         north_grid[min_inlat+2][min_inlon], 
                         current_easting, current_northing, &move_down);
                 else
                     move_down = (float) SHRT_MAX;

                 if ((move_right - stay_right) < MINSIGMA)
                    min_inlon++;
                 if ((move_up - stay_up) < MINSIGMA)
                    min_inlat--;
                 if ((move_down - stay_down) < MINSIGMA)
                    min_inlat++;

                    /* extract UTM coordinates of four points to be interpolated 
                       and build array */
                    coordinates[0][0] =  east_grid[min_inlat][min_inlon];
                    coordinates[0][1] =  north_grid[min_inlat][min_inlon];
                    coordinates[1][0] =  east_grid[min_inlat][min_inlon+1];
                    coordinates[1][1] =  north_grid[min_inlat][min_inlon+1];
                    coordinates[2][0] =  east_grid[min_inlat+1][min_inlon];
                    coordinates[2][1] =  north_grid[min_inlat+1][min_inlon];
                    coordinates[3][0] =  east_grid[min_inlat+1][min_inlon+1];
                    coordinates[3][1] =  north_grid[min_inlat+1][min_inlon+1];

                    /* determine index of four points in order to pull from 
                       atmospheric parameter file */
                    indices[0] = min_inlat * num_jays + min_inlon;
                    indices[1] = (min_inlat + 1) * num_jays + min_inlon;
                    indices[2] = min_inlat * num_jays + min_inlon + 1;
                    indices[3] = (min_inlat + 1) * num_jays + min_inlon + 1;
                }

                /* convert height from m to km */
                dem[col] = (float) dem[col] / 1000.0;

                /* interpolate three parameters to that height at each of the 
                   four closest points */
                for (g = 0; g < 4; g++)
                {
                    current_index = indices[g] * NUM_ELEVATIONS;
                    /* extract atmospheric parameters for all heights at current 
                       location*/
                    for (n = 0; n < 6; n++)
                    {
                        counter = 0;
                        for (k = current_index; k < current_index + 
                            NUM_ELEVATIONS - 1; k++)
                        {
                                current_location[n][counter] = 
                                    modtran_results[n][k];
                        }
                        counter++;
                    }

                    /* separate height and atmospheric data */
                    for (k = current_index; k < current_index + 4; k++)
                    {
                        /* interpolate three atmospheric parameters to current 
                           height*/
                        interpolate_to_height(k, current_location[2], 
                            current_location[3], current_location[4], 
                            current_location[5], dem[col], at_height[g]);
                    }
                }

                /* interpolate parameters at appropriate height to location of 
                   current pixel*/
                interpolate_to_location(coordinates, at_height, current_easting,
                                        current_northing, parameters);

                /* convert radiances to W*m^(-2)*sr(-1) */
                landsat_results[0][col] = parameters[0];
                landsat_results[1][col] = 10000.0 * parameters[1];
                landsat_results[2][col] = 10000.0 * parameters[2];

            }
        }

        /* Write out the temporary binary output files
           Note: needs to be deleted before release */
        status = fwrite(&input->therm_buf[0], sizeof(int16), 
                 input->size_th.s, therm_fptr);
        if (status != input->size_th.s)
        {
            sprintf(errstr, "Writing to %s", therm_fname);
            LST_ERROR (errstr, "third_pixels_post");
        }

        status = fwrite(&landsat_results[0][0], sizeof(float), 
                 input->size_th.s, trans_fptr);
        if (status != input->size_th.s)
        {
            sprintf(errstr, "Writing to %s", trans_fname);
            LST_ERROR (errstr, "third_pixels_post");
        }

        status = fwrite(&landsat_results[1][0], sizeof(float), 
                 input->size_th.s, up_fptr);
        if (status != input->size_th.s)
        {
            sprintf(errstr, "Writing to %s", up_fname);
            LST_ERROR (errstr, "third_pixels_post");
        }

        status = fwrite(&landsat_results[2][0], sizeof(float), 
                 input->size_th.s, down_fptr);
        if (status != input->size_th.s)
        {
            sprintf(errstr, "Writing to %s", down_fname);
            LST_ERROR (errstr, "third_pixels_post");
        }
    }

    /* Free allocated memory */
    free(narr_lat);
    free(narr_lon);
    free(distances);
    free(distance_copy);
    status = free_2d_array((void **)at_height);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: at_height\n");
        LST_ERROR (errstr, "third_pixels_post");              
    }
    status = free_2d_array((void **)coordinates);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: coordinates\n");
        LST_ERROR (errstr, "third_pixels_post");              
    }
    status = free_2d_array((void **)east_grid);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: east_grid\n");
        LST_ERROR (errstr, "third_pixels_post");              
    }
    status = free_2d_array((void **)north_grid);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: north_grid\n");
        LST_ERROR (errstr, "third_pixels_post");              
    }

    status = free_2d_array((void **)narr_utm);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: narr_utm\n");
        LST_ERROR (errstr, "third_pixels_post");              
    }

    status = free_2d_array((void **)landsat_results);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: landsat_results\n");
        LST_ERROR (errstr, "third_pixels_post");              
    }

    /* Close the intermediate binary files
       Note: needs to be deleted before release */
    status = fclose(therm_fptr);
    if ( status )
    {
        sprintf(errstr, "Closing file %s", therm_fname);
        LST_ERROR (errstr, "third_pixels_post");
    }

    status = fclose(trans_fptr);
    if ( status )
    {
        sprintf(errstr, "Closing file %s", trans_fname);
        LST_ERROR (errstr, "third_pixels_post");
    }

    status = fclose(up_fptr);
    if ( status )
    {
        sprintf(errstr, "Closing file %s", up_fname);
        LST_ERROR (errstr, "third_pixels_post");
    }

    status = fclose(down_fptr);
    if ( status )
    {
        sprintf(errstr, "Closing file %s", down_fname);
        LST_ERROR (errstr, "third_pixels_post");
    }

    return SUCCESS;
}
