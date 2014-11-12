#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <sys/stat.h>
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

/******************************************************************************
MODULE:  convert_geopotential_geometric

PURPOSE: Convert array of geopotential heights to array of geometric heights 
         given latitude 

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/22/2014   Song Guo         Original Development
******************************************************************************/
int convert_geopotential_geometric
(
    int num_points,
    float *lat,
    float **geo_potential,
    float **geo_metric
)
{
    float *radlat;
    int i, j;
    int num_pressures = P_LAYER;
    float g_0 = 9.80665;
    float r_max = 6378.137;
    float r_min = 6356.752;
    float *radius;
    float *gravity_ratio;
    char errstr[MAX_STR_LEN];

    /* Allocate memeory*/
    radlat = (float *)malloc(num_points * sizeof(float)); 
    if (radlat == NULL)
    {
        sprintf (errstr, "Allocating radlat memory");
        LST_ERROR (errstr, "convert_geopotential_geometric");
    }

    radius = (float *)malloc(num_points * sizeof(float)); 
    if (radius == NULL)
    {
        sprintf (errstr, "Allocating radius memory");
        LST_ERROR (errstr, "convert_geopotential_geometric");
    }

    gravity_ratio = (float *)malloc(num_points * sizeof(float)); 
    if (gravity_ratio == NULL)
    {
        sprintf (errstr, "Allocating gravity_ratio memory");
        LST_ERROR (errstr, "convert_geopotential_geometric");
    }

    for (i = 0; i < num_points; i++)
    {
        radlat[i] = (lat[i] * PI) / 180.0;

        /* define variable based on latitude */
        radius[i] = 1000.0 * (sqrt(1.0/((cos(radlat[i])*cos(radlat[i]))/(r_max*r_max)+
                    ((sin(radlat[i])*sin(radlat[i]))/(r_min*r_min)))));
        gravity_ratio[i] = (9.80616*(1-0.002637*cos(2*radlat[i])+0.0000059*
               (cos(2*radlat[i])*cos(2*radlat[i]))))/g_0;
    }

    /* Dynamic allocate the geo_metric memory */
    geo_metric = (float **)allocate_2d_array(num_pressures, NARR_ROW * NARR_COL,  
                 sizeof(float)); 
    if (geo_metric == NULL)
    {
        sprintf (errstr, "Allocating geo_metric memory");
        LST_ERROR (errstr, "convert_geopotential_geometric");
    }

    for (i = 0; i < P_LAYER; i++)
    {
        for ( j = 0; j < num_points; j++)
        {
            geo_metric[i][j] = (geo_potential[i][j] * radius[j]) / (1000.0 * 
                  (gravity_ratio[j] *radius[j] - geo_potential[i][j]));
        }
    }

    /* free memory */
    free(radlat);
    free(radius);
    free(gravity_ratio); 

    return SUCCESS;
}
/******************************************************************************
MODULE:  convert_sh_rh

PURPOSE: Given array of specific humidities, temperature, and pressure, generate 
         array of relative humidities

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/22/2014   Song Guo         Original Development
******************************************************************************/
int convert_sh_rh
(
    int num_points,
    float *lat,
    float **spec_hum,
    float **temp_k,
    float **pressure,
    float **rh
)
{
    char errstr[MAX_STR_LEN];
    int i, j;
    float mh20 = 18.01534;
    float mdry = 28.9644;

    float a0w = 6.107799961;
    float a1w = 4.436518521e-1;
    float a2w = 1.428945805e-2;
    float a3w = 2.650648471e-4;
    float a4w = 3.0312403963e-6;
    float a5w = 2.034080948e-8;
    float a6w = 6.136820929e-11;
    float **temp_c;
    float **ewater;
    float **e2;
    float **goff;
    float **ph20;
    int status;

    /* Allocate memory */
    temp_c = (float **)allocate_2d_array(P_LAYER, num_points, sizeof(float)); 
    if (temp_c == NULL)
    {
        sprintf (errstr, "Allocating temp_c memory");
        LST_ERROR (errstr, "first_files");
    }

    ewater = (float **)allocate_2d_array(P_LAYER, num_points, sizeof(float)); 
    if (ewater == NULL)
    {
        sprintf (errstr, "Allocating ewater memory");
        LST_ERROR (errstr, "first_files");
    }

    e2 = (float **)allocate_2d_array(P_LAYER, num_points, sizeof(float)); 
    if (e2 == NULL)
    {
        sprintf (errstr, "Allocating e2 memory");
        LST_ERROR (errstr, "first_files");
    }

    goff = (float **)allocate_2d_array(P_LAYER, num_points, sizeof(float)); 
    if (goff == NULL)
    {
        sprintf (errstr, "Allocating goff memory");
        LST_ERROR (errstr, "first_files");
    }

    ph20 = (float **)allocate_2d_array(P_LAYER, num_points, sizeof(float)); 
    if (ph20 == NULL)
    {
        sprintf (errstr, "Allocating  memory");
        LST_ERROR (errstr, "first_files");
    }


    for (i = 0; i < P_LAYER; i++)
    {
        for ( j = 0; j < num_points; j++)
        {
            /* Convert temperature to C */
            temp_c[i][j] = temp_k[i][j] - 273.15;
            /* calculate vapor pressure at given temperature */
            ewater[i][j] = a0w + temp_c[i][j] * (a1w + temp_c[i][j] * (a2w + 
                temp_c[i][j] * (a3w + temp_c[i][j] * (a4w + temp_c[i][j] * 
                (a5w + temp_c[i][j] * (a6w * temp_c[i][j])))))); /* hpa */

            e2[i][j] = exp(-0.58002206e4/temp_k[i][j] + 0.13914993 - 0.48640239e-1 * 
                 temp_k[i][j]+0.41764768e-4 * pow(temp_k[i][j], 2) - 0.14452093e-7 * 
                 pow(temp_k[i][j], 3) + 0.65459673 * log(temp_k[i][j])); /* Pa */

            goff[i][j] = -7.90298 * (373.16/temp_k[i][j]-1) + 5.02808 * log10(373.16 / 
                temp_k[i][j]) - 1.3816e-7 * pow(10, (11.344 * (1 - (temp_k[i][j] / 373.16)))
                - 1) + 8.1328e-3 * pow(10, (-3.49149 * (373.16 / temp_k[i][j] - 1)) - 1)
                + log10(1013.246); /* hPa */

             ph20[i][j] = (spec_hum[i][j] * pressure[i][j] * mdry) / (mh20 - 
                 spec_hum[i][j] * mh20 + spec_hum[i][j] * mdry); 

             rh[i][j] = (ph20[i][j] / pow(10, goff[i][j])) * 100.0;
        }
    } 

    /* Free allocated memory */
    status = free_2d_array((void **)temp_c);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: temp_c\n");
        LST_ERROR (errstr, "first_files");              
    }

    status = free_2d_array((void **)ewater);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: ewater\n");
        LST_ERROR (errstr, "first_files");              
    }

    status = free_2d_array((void **)e2);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: e2\n");
        LST_ERROR (errstr, "first_files");              
    }

    status = free_2d_array((void **)goff);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: goff\n");
        LST_ERROR (errstr, "first_files");              
    }

    status = free_2d_array((void **)ph20);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: ph20\n");
        LST_ERROR (errstr, "first_files");              
    }

    return SUCCESS;
}

/******************************************************************************
MODULE:  first_files

PURPOSE: Creates directories and writes tape5 file, caseList, and commandList

RETURN: SUCCESS
        FAILURE

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
9/21/2014   Song Guo         Original Development
******************************************************************************/
int first_files
(
    Input_t *input,             /*I: input structure */
    char **case_list,           /*O: modtran run list */
    char **command_list,        /*O: modtran run command list */
    int *entry,                 /*O: number of cases/commands */
    int *num_points,            /*O: number of NARR points */
    bool verbose                /*I: value to indicate if intermediate messages 
                                     be printed */
)
{
    char errstr[MAX_STR_LEN];  /* error string */
    int **eye;
    int **jay;
    float **lat;
    float **lon;
    float **hgt1;
    float **shum1;
    float **tmp1;
    float **hgt2;
    float **shum2;
    float **tmp2;
    float *narr_lat;
    float *narr_lon;
    float **narr_hgt1;
    float **narr_shum1;
    float **narr_tmp1;
    float **narr_hgt2;
    float **narr_shum2;
    float **narr_tmp2;
    float **narr_height;
    float **narr_height1;
    float **narr_height2;
    float **narr_rh;
    float **narr_rh1;
    float **narr_rh2;
    float **narr_tmp;
    int i, j, k;
    int p[P_LAYER] = {1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 725, 700, 
        650, 600, 550, 500, 450, 400, 350, 300, 275, 250, 225, 200, 175, 150, 125, 100};
    char full_path[MAX_STR_LEN];
    int landsat_hemi;
    float narr_ul_lat;
    float narr_ul_lon;
    float narr_lr_lat;
    float narr_lr_lon;
    int in_counter = 0;
    int max_eye;
    int min_eye;
    int max_jay;
    int min_jay;
    int num_eyes;
    int num_jays;
    float **pressure;
    int rem1;
    int rem2;
    float hour1;
    float hour2;
    float time;
    FILE *fd;
    float *stan_height;
    float *stan_pre;
    float *stan_temp;
    float *stan_rh;
    float *temp_height;
    float *temp_pressure;
    float *temp_temp;
    float *temp_rh;
    float gndalt[NUM_ELEVATIONS] = {0.0, 0.6, 1.1, 1.6, 2.1, 2.6, 3.1, 3.6, 4.05};
    int num_cases;
    char command[MAX_STR_LEN];
    char current_gdalt[MAX_STR_LEN];
    char current_temp[MAX_STR_LEN];
    char current_alb[MAX_STR_LEN];
    char current_point[MAX_STR_LEN];
    char temp_out[MAX_STR_LEN];
    int index_below = 0;
    int index_above = NUM_ELEVATIONS;
    float new_height;
    float new_pressure;
    float new_temp;
    float new_rh;
    int index, index2;
    int *counter;
    float tmp[3] = {273.0, 310.0, 0.0};
    float alb[3] = {0.0, 0.0, 0.1};
    char *path = NULL; 
    char *modtran_path = NULL; 
    struct stat s;
    int status;
    int temp_int1, temp_int2; 

    /* Dynamic allocate the 2d memory */
    eye = (int **)allocate_2d_array(NARR_ROW, NARR_COL, sizeof(int)); 
    if (eye == NULL)
    {
        sprintf (errstr, "Allocating eye memory");
        LST_ERROR (errstr, "first_files");
    }

    jay = (int **)allocate_2d_array(NARR_ROW, NARR_COL, sizeof(int)); 
    if (jay == NULL)
    {
        sprintf (errstr, "Allocating jay memory");
        LST_ERROR (errstr, "first_files");
    }

    lat = (float **)allocate_2d_array(NARR_ROW, NARR_COL, sizeof(float)); 
    if (lat == NULL)
    {
        sprintf (errstr, "Allocating lat memory");
        LST_ERROR (errstr, "first_files");
    }

    lon = (float **)allocate_2d_array(NARR_ROW, NARR_COL, sizeof(float)); 
    if (lon == NULL)
    {
        sprintf (errstr, "Allocating lon memory");
        LST_ERROR (errstr, "first_files");
    }

    path = getenv("LST_DATA");
    if (path == NULL)
    {
        sprintf (errstr, "LST_DATA environment variable is not set");
        LST_ERROR(errstr, "first_files");
    }

    sprintf(full_path,"%s/%s",path,"coordinates.txt");
    fd = fopen(full_path, "r");
    if (fd == NULL)
    {
        sprintf (errstr, "Can't open coordinates.txt file");
        LST_ERROR(errstr, "first_files");
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
            if ((lon[i][j] - 180.0) > MINSIGMA)
                lon[i][j] = 360.0 - lon[i][j];
            else
                lon[i][j] = - lon[i][j];
        }
    }
    fclose(fd);

    /* Dynamic allocate the 2d memory */
    hgt1 = (float **)allocate_2d_array(P_LAYER, NARR_ROW * NARR_COL, sizeof(float)); 
    if (hgt1 == NULL)
    {
        sprintf (errstr, "Allocating hgt_1 memory");
        LST_ERROR (errstr, "first_files");
    }

    shum1 = (float **)allocate_2d_array(P_LAYER, NARR_ROW * NARR_COL, sizeof(float)); 
    if (shum1 == NULL)
    {
        sprintf (errstr, "Allocating shum_1 memory");
        LST_ERROR (errstr, "first_files");
    }

    tmp1 = (float **)allocate_2d_array(P_LAYER, NARR_ROW * NARR_COL, sizeof(float)); 
    if (tmp1 == NULL)
    {
        sprintf (errstr, "Allocating tmp_1 memory");
        LST_ERROR (errstr, "first_files");
    }

    /* Read in NARR height for time before landsat acqusition */ 
    for (i = 0; i < P_LAYER - 1; i++)
    {
        sprintf(full_path,"%s/%d%s","HGT_1/", p[i], ".txt");
        fd = fopen(full_path, "r");
        if (fd == NULL)
        {
            sprintf (errstr, "Can't HGT_1 txt file");
            LST_ERROR(errstr, "first_files");
        }

        fscanf(fd, "%d %d", &temp_int1, &temp_int2);
        for (j = 0; j < NARR_ROW * NARR_COL - 1; j++)
        {
            if (fscanf(fd, "%f", &hgt1[i][j]) == EOF)
            {
                sprintf (errstr, "End of file (EOF) is met before "
                         "NARR_ROW * NARR_COL lines");
                LST_ERROR(errstr, "first_files");
            }
        }
        fclose(fd);
    }

    /* Read in NARR specific humidity for time before landsat acqusition */ 
    for (i = 0; i < P_LAYER - 1; i++)
    {
        sprintf(full_path,"%s/%d%s","SHUM_1/", p[i], ".txt");
        fd = fopen(full_path, "r");
        if (fd == NULL)
        {
            sprintf (errstr, "Can't open SHUM_1 file");
            LST_ERROR(errstr, "first_files");
        }

        fscanf(fd, "%d %d", &temp_int1, &temp_int2);
        for (j = 0; j < NARR_ROW * NARR_COL - 1; j++)
        {
            if (fscanf(fd, "%f", &shum1[i][j]) == EOF)
            {
                sprintf (errstr, "End of file (EOF) is met before "
                         "NARR_ROW * NARR_COL lines");
                LST_ERROR(errstr, "first_files");
            }
        }
        fclose(fd);
    }

    /* Read in NARR temperature for time before landsat acqusition */ 
    for (i = 0; i < P_LAYER - 1; i++)
    {
        sprintf(full_path,"%s/%d%s","TMP_1/", p[i], ".txt");
        fd = fopen(full_path, "r");
        if (fd == NULL)
        {
            sprintf (errstr, "Can't open TMP_1 file");
            LST_ERROR(errstr, "first_files");
        }

        fscanf(fd, "%d %d", &temp_int1, &temp_int2);
        for (j = 0; j < NARR_ROW * NARR_COL - 1; j++)
        {
            if (fscanf(fd, "%f", &tmp1[i][j]) == EOF)
            {
                sprintf (errstr, "End of file (EOF) is met before "
                         "NARR_ROW * NARR_COL lines");
                LST_ERROR(errstr, "first_files");
            }
        }
        fclose(fd); 
    }

    /* Dynamic allocate the 2d memory */
    hgt2 = (float **)allocate_2d_array(P_LAYER, NARR_ROW * NARR_COL, sizeof(float)); 
    if (hgt2 == NULL)
    {
        sprintf (errstr, "Allocating hgt_2 memory");
        LST_ERROR (errstr, "first_files");
    }

    shum2 = (float **)allocate_2d_array(P_LAYER, NARR_ROW * NARR_COL, sizeof(float)); 
    if (shum2 == NULL)
    {
        sprintf (errstr, "Allocating shum_2 memory");
        LST_ERROR (errstr, "first_files");
    }

    tmp2 = (float **)allocate_2d_array(P_LAYER, NARR_ROW * NARR_COL, sizeof(float)); 
    if (tmp2 == NULL)
    {
        sprintf (errstr, "Allocating tmp_2 memory");
        LST_ERROR (errstr, "first_files");
    }

    /* Read in NARR height for time after landsat acqusition */ 
    for (i = 0; i < P_LAYER - 1; i++)
    {
        sprintf(full_path,"%s/%d%s","HGT_2/", p[i], ".txt");
        fd = fopen(full_path, "r");
        if (fd == NULL)
        {
            sprintf (errstr, "Can't HGT_2 txt file");
            LST_ERROR(errstr, "first_files");
        }

        fscanf(fd, "%d %d", &temp_int1, &temp_int2);
        for (j = 0; j < NARR_ROW * NARR_COL - 1; j++)
        {
            if (fscanf(fd, "%f", &hgt2[i][j]) == EOF)
            {
                sprintf (errstr, "End of file (EOF) is met before "
                         "NARR_ROW * NARR_COL lines");
                LST_ERROR(errstr, "first_files");
            }
        }
        fclose(fd);
    }
 
    /* Read in NARR specific humidity for time after landsat acqusition */ 
    for (i = 0; i < P_LAYER - 1; i++)
    {
        sprintf(full_path,"%s/%d%s","SHUM_2/", p[i], ".txt");
        fd = fopen(full_path, "r");
        if (fd == NULL)
        {
            sprintf (errstr, "Can't open SHUM_2 file");
            LST_ERROR(errstr, "first_files");
        }

        fscanf(fd, "%d %d", &temp_int1, &temp_int2);
        for (j = 0; j < NARR_ROW * NARR_COL - 1; j++)
        {
            if (fscanf(fd, "%f", &shum2[i][j]) == EOF)
            {
                sprintf (errstr, "End of file (EOF) is met before "
                         "NARR_ROW * NARR_COL lines");
                LST_ERROR(errstr, "first_files");
            }
        }
        fclose(fd);
    }

    /* Read in NARR temperature for time after landsat acqusition */ 
    for (i = 0; i < P_LAYER - 1; i++)
    {
        sprintf(full_path,"%s/%d%s","TMP_2/", p[i], ".txt");
        fd = fopen(full_path, "r");
        if (fd == NULL)
        {
            sprintf (errstr, "Can't open TMP_2 file");
            LST_ERROR(errstr, "first_files");
        }

        fscanf(fd, "%d %d", &temp_int1, &temp_int2);
        for (j = 0; j < NARR_ROW * NARR_COL - 1; j++)
        {
            if (fscanf(fd, "%f", &tmp2[i][j]) == EOF)
            {
                sprintf (errstr, "End of file (EOF) is met before "
                         "NARR_ROW * NARR_COL lines");
                LST_ERROR(errstr, "first_files");
            }
        }
        fclose(fd); 
    }

    /* determine if landsat is in the northern or southern hemisphere. 
       '6' = northern hemisphere, '7' = southern hermisphere. */
    if (input->meta.ul_geo_corner.lat > MINSIGMA)
        landsat_hemi = 6;
    else
        landsat_hemi = 7; 

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
    *num_points = num_eyes * num_jays;

    /* Allocate memory for height of NARR points within the rectangular */
    narr_lat = (float *)malloc(*num_points * sizeof(float)); 
    if (narr_lat == NULL)
    {
        sprintf (errstr, "Allocating narr_lat memory");
        LST_ERROR (errstr, "first_files");
    }

    narr_lon = (float *)malloc(*num_points * sizeof(float)); 
    if (narr_lon == NULL)
    {
        sprintf (errstr, "Allocating narr_lon memory");
        LST_ERROR (errstr, "first_files");
    }

    narr_hgt1 = (float **)allocate_2d_array(P_LAYER, *num_points,sizeof(float)); 
    if (narr_hgt1 == NULL)
    {
        sprintf (errstr, "Allocating narr_hgt1 memory");
        LST_ERROR (errstr, "first_files");
    }

    narr_hgt2 = (float **)allocate_2d_array(P_LAYER, *num_points, sizeof(float)); 
    if (narr_hgt2 == NULL)
    {
        sprintf (errstr, "Allocating narr_hgt2 memory");
        LST_ERROR (errstr, "first_files");
    }

    /* Allocate memory for humidity of NARR points within the rectangular */
    narr_shum1 = (float **)allocate_2d_array(P_LAYER, *num_points, sizeof(float)); 
    if (narr_shum1 == NULL)
    {
        sprintf (errstr, "Allocating narr_shum1 memory");
        LST_ERROR (errstr, "first_files");
    }

    narr_shum2 = (float **)allocate_2d_array(P_LAYER, *num_points, sizeof(float)); 
    if (narr_shum2 == NULL)
    {
        sprintf (errstr, "Allocating narr_shum2 memory");
        LST_ERROR (errstr, "first_files");
    }

    /* Allocate memory for temperature of NARR points within the rectangular */
    narr_tmp1 = (float **)allocate_2d_array(P_LAYER, *num_points, sizeof(float)); 
    if (narr_tmp1 == NULL)
    {
        sprintf (errstr, "Allocating narr_tmp1 memory");
        LST_ERROR (errstr, "first_files");
    }

    narr_tmp2 = (float **)allocate_2d_array(P_LAYER, *num_points, sizeof(float)); 
    if (narr_tmp2 == NULL)
    {
        sprintf (errstr, "Allocating narr_tmp2 memory");
        LST_ERROR (errstr, "first_files");
    }

    /* extract coordinates within the NARR rectangle */
    for (i = min_eye; i <= max_eye; i++)
    {
        for (j = min_jay; j <= max_jay; j++)
        {
            narr_lat[(i-min_eye) * num_jays + (j-min_jay)] = lat[i][j]; 
            narr_lon[(i-min_eye) * num_jays + (j-min_jay)] = lon[i][j]; 
        }
    }
    for (k = 0; k < P_LAYER; k++)
    {
        for (i = min_eye; i <= max_eye; i++)
        {
            for (j = min_jay; j <= max_jay; j++)
            {
                narr_hgt1[k][(i-min_eye) * num_jays + (j-min_jay)] 
                      = hgt1[k][i * num_jays + j]; 
                narr_shum1[k][(i-min_eye) * num_jays + (j-min_jay)]  
                      = shum1[k][i * num_jays + j]; 
                narr_tmp1[k][(i-min_eye) * num_jays + (j-min_jay)]  
                      = tmp1[k][i * num_jays + j]; 
                narr_hgt2[k][(i-min_eye) * num_jays + (j-min_jay)]  
                      = hgt2[k][i * num_jays + j]; 
                narr_shum2[k][(i-min_eye) * num_jays + (j-min_jay)]  
                      = shum2[k][i * num_jays + j]; 
                narr_tmp2[k][(i-min_eye) * num_jays + (j-min_jay)] 
                      = tmp2[k][i * num_jays + j]; 
            }
        }
    }

    /* Release memory */
    status = free_2d_array((void **)eye);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: eye\n");
        LST_ERROR (errstr, "first_files");              
    }

    status = free_2d_array((void **)jay);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: jay\n");
        LST_ERROR (errstr, "first_files");              
    }

    status = free_2d_array((void **)lat);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: lat\n");
        LST_ERROR (errstr, "first_files");              
    }

    status = free_2d_array((void **)lon);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: lon\n");
        LST_ERROR (errstr, "first_files");              
    }

    status = free_2d_array((void **)hgt1);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: hgt1\n");
        LST_ERROR (errstr, "first_files");              
    }

    status = free_2d_array((void **)shum1);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: shum1\n");
        LST_ERROR (errstr, "first_files");              
    }

    status = free_2d_array((void **)tmp1);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: tmp1\n");
        LST_ERROR (errstr, "first_files");              
    }

    status = free_2d_array((void **)hgt2);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: ght2\n");
        LST_ERROR (errstr, "first_files");              
    }

    status = free_2d_array((void **)shum2);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: shum2\n");
        LST_ERROR (errstr, "first_files");              
    }

    status = free_2d_array((void **)tmp2);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: tmp2\n");
        LST_ERROR (errstr, "first_files");              
    }

    /* Allocate memory */
    pressure = (float **)allocate_2d_array(P_LAYER, *num_points, sizeof(float)); 
    if (pressure == NULL)
    {
        sprintf (errstr, "Allocating pressure memory");
        LST_ERROR (errstr, "first_files");
    }

    narr_height1 = (float **)allocate_2d_array(P_LAYER, *num_points, sizeof(float)); 
    if (narr_height1 == NULL)
    {
        sprintf (errstr, "Allocating narr_height1 memory");
        LST_ERROR (errstr, "first_files");
    }

    narr_height2 = (float **)allocate_2d_array(P_LAYER, *num_points, sizeof(float)); 
    if (narr_height2 == NULL)
    {
        sprintf (errstr, "Allocating narr_height2 memory");
        LST_ERROR (errstr, "first_files");
    }

    narr_rh1 = (float **)allocate_2d_array(P_LAYER, *num_points, sizeof(float)); 
    if (narr_rh1 == NULL)
    {
        sprintf (errstr, "Allocating narr_rh1 memory");
        LST_ERROR (errstr, "first_files");
    }

    narr_rh2 = (float **)allocate_2d_array(P_LAYER, *num_points, sizeof(float)); 
    if (narr_rh2 == NULL)
    {
        sprintf (errstr, "Allocating narr_rh2 memory");
        LST_ERROR (errstr, "first_files");
    }

    for (i = 0; i < P_LAYER; i++)
    {
        for ( j = 0; j < *num_points; j++)
        {
            pressure[i][j] = p[i];
        }
    }

    /* convert grib data to variables to be input to MODTRAN */
    status = convert_geopotential_geometric(*num_points, narr_lat, narr_hgt1, 
                                            narr_hgt1);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Calling convert_geopotential_geometric1");
        LST_ERROR (errstr, "first_files");
    }
 
    status = convert_geopotential_geometric(*num_points, narr_lat, narr_hgt2, 
                                            narr_hgt2);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Calling convert_geopotential_geometric2");
        LST_ERROR (errstr, "first_files");
    }
 
    status = convert_sh_rh(*num_points, narr_lat, narr_shum1, narr_tmp1, 
                           pressure, narr_rh1);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Calling convert_sh_rh1");
        LST_ERROR (errstr, "first_files");
    }
 
    status = convert_sh_rh(*num_points, narr_lat, narr_shum2, narr_tmp2, 
                           pressure, narr_rh2);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Calling convert_sh_rh1");
        LST_ERROR (errstr, "first_files");
    }

    /* free allocated memory */
    status = free_2d_array((void **)narr_hgt1);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: narr_hgt1\n");
        LST_ERROR (errstr, "first_files");              
    }

    status = free_2d_array((void **)narr_hgt2);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: narr_hgt2\n");
        LST_ERROR (errstr, "first_files");              
    }

    status = free_2d_array((void **)narr_shum1);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: narr_shum1\n");
        LST_ERROR (errstr, "first_files");              
    }

    status = free_2d_array((void **)narr_shum2);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: narr_shum2\n");
        LST_ERROR (errstr, "first_files");              
    }

    /* determine three hour-increment before and after scene center scan time */
    rem1 = input->meta.acq_date.hour % 3;
    rem2 = 3 - rem1;
    hour1 = (float) (input->meta.acq_date.hour - rem1);
    hour2 = (float) (input->meta.acq_date.hour + rem2);

    /* Round to the nearest minute */
    if ((input->meta.acq_date.second -30.0) >= MINSIGMA)
        input->meta.acq_date.minute++; 
 
    /* convert hour-min acquisition time to decimal time */
    time = (float)input->meta.acq_date.hour + (float) input->meta.acq_date.minute 
           / 60.0;

    /* Allocate memory */
    narr_height = (float **)allocate_2d_array(P_LAYER, *num_points, sizeof(float)); 
    if (narr_height == NULL)
    {
        sprintf (errstr, "Allocating narr_height memory");
        LST_ERROR (errstr, "first_files");
    }

    /* Allocate memory */
    narr_rh = (float **)allocate_2d_array(P_LAYER,*num_points, sizeof(float)); 
    if (narr_rh == NULL)
    {
        sprintf (errstr, "Allocating narr_rh memory");
        LST_ERROR (errstr, "first_files");
    }

    /* Allocate memory */
    narr_tmp = (float **)allocate_2d_array(P_LAYER, *num_points, sizeof(float)); 
    if (narr_tmp == NULL)
    {
        sprintf (errstr, "Allocating narr_tmp memory");
        LST_ERROR (errstr, "first_files");
    }

    /* linearly interpolate geometric height, relative humidity, and temperature 
       for NARR points within Landsat scene this is the NARR data corresponding 
       to the acquisition time of the Landsat image converted to appropriated
       variable for MODTRAN input */
    for (i = 0; i < P_LAYER; i++)
    {
        for ( j = 0; j < *num_points; j++)
        {
            narr_height[i][j] = narr_height1[i][j] + (time - hour1) * 
                ((narr_height2[i][j] - narr_height1[i][j]) / (hour2 - hour1)); 
            narr_rh[i][j] = narr_rh1[i][j] + (time - hour1) * 
                ((narr_rh2[i][j] - narr_rh1[i][j]) / (hour2 - hour1)); 
            narr_tmp[i][j] = narr_tmp1[i][j] + (time - hour1) * 
                ((narr_tmp2[i][j] - narr_tmp1[i][j]) / (hour2 - hour1)); 
        }
    }

    /* Free allocated memory */
    status = free_2d_array((void **)narr_height1);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: narr_height1\n");
        LST_ERROR (errstr, "first_files");              
    }
    
    status = free_2d_array((void **)narr_height2);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: narr_height2\n");
        LST_ERROR (errstr, "first_files");              
    }
    
    status = free_2d_array((void **)narr_rh1);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: narr_rh1\n");
        LST_ERROR (errstr, "first_files");              
    }
    
    status = free_2d_array((void **)narr_rh2);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: narr_rh2\n");
        LST_ERROR (errstr, "first_files");              
    }
    
    status = free_2d_array((void **)narr_tmp1);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: narr_tmp1\n");
        LST_ERROR (errstr, "first_files");              
    }
    
    status = free_2d_array((void **)narr_tmp2);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: narr_tmp2\n");
        LST_ERROR (errstr, "first_files");              
    }

    /* Build tape 5 files */
    /* Allocate memory */
    stan_height = (float *)malloc(STAN_LAYER * sizeof(float)); 
    if (stan_height == NULL)
    {
        sprintf (errstr, "Allocating stan_height memory");
        LST_ERROR (errstr, "first_files");
    }

    stan_pre = (float *)malloc(STAN_LAYER * sizeof(float)); 
    if (stan_pre == NULL)
    {
        sprintf (errstr, "Allocating stan_pre memory");
        LST_ERROR (errstr, "first_files");
    }

    stan_temp = (float *)malloc(STAN_LAYER * sizeof(float)); 
    if (stan_temp == NULL)
    {
        sprintf (errstr, "Allocating stan_temp memory");
        LST_ERROR (errstr, "first_files");
    }

    stan_rh = (float *)malloc(STAN_LAYER * sizeof(float)); 
    if (stan_rh == NULL)
    {
        sprintf (errstr, "Allocating stan_rh memory");
        LST_ERROR (errstr, "first_files");
    }

    /* read in file containing standard mid lat summer atmosphere information to be 
       used for upper layers */
    sprintf(full_path,"%s/%s",path,"stanAtm.txt");
    fd = fopen(full_path, "r"); 
    if (fd == NULL)
    {
        sprintf (errstr, "Opening file: stanAtm.txt\n");
        LST_ERROR (errstr, "first_files");
    }

    for (i = 0; i < STAN_LAYER; i++)
    {
        if (fscanf(fd, "%f,%f,%f,%f", &stan_height[i], &stan_pre[i], &stan_temp[i], 
                &stan_rh[i]) == EOF)
        {
            sprintf (errstr, "End of file (EOF) is met before STAN_LAYER lines");
            LST_ERROR(errstr, "first_files");
        }
    }

    status = fclose(fd);
    if ( status )
    {
        sprintf (errstr, "Closing file: stanAtm.txt\n");
        LST_ERROR (errstr, "first_files");
    }

    /* determine number of MODTRAN runs */
    num_cases = *num_points * NUM_ELEVATIONS * 3;  

    /* Allocate memory */
    counter = (int *)malloc(STAN_LAYER * sizeof(int));
    if (counter == NULL)
    {
        sprintf (errstr, "Allocating counter memory");
        LST_ERROR (errstr, "first_files");
    }

    case_list = (char **)allocate_2d_array(num_cases, MAX_STR_LEN, sizeof(char));  
    if (case_list == NULL)
    {
        sprintf (errstr, "Allocating case_list memory");
        LST_ERROR (errstr, "first_files");
    }

    command_list = (char **)allocate_2d_array(num_cases, MAX_STR_LEN, sizeof(char));  
    if (command_list == NULL)
    {
        sprintf (errstr, "Allocating command_list memory");
        LST_ERROR (errstr, "first_files");
    }

    modtran_path = getenv("MODTRAN_PATH");
    if (modtran_path == NULL)
    {
        sprintf (errstr, "MODTRAN_PATH environment variable is not set");
        LST_ERROR(errstr, "first_files");
    }

    for (i = 0; i < *num_points; i++)
    {
        /* create a directory for the current NARR point */
        if (narr_lon[i] < 0)
            narr_lon[i] = -narr_lon[i];
        else
            narr_lon[i] = 360.0 - narr_lon[i];

        if ((narr_lon[i] - 100.0) < MINSIGMA)
            sprintf(current_point,"%6.3f_%6.3f", narr_lat[i], narr_lon[i]);
        else
            sprintf(current_point,"%6.3f_%6.2f", narr_lat[i], narr_lon[i]); 

        if (stat(current_point, &s)== -1)
        {
            status = mkdir(current_point, 0755);
            if (status != SUCCESS)
            {
                sprintf (errstr, "system call 1");
                LST_ERROR (errstr, "first_files");
            }
        }

        /* set lowest altitude is the first geometric height at that NARR point 
           (if positive) and (if negative set to zero) */
        if (narr_height[0][i] < 0)
            gndalt[0] = 0.0;   
        else
            gndalt[0] = narr_height[0][i];

        /* determine latitude and longitude of current NARR point and insert into 
           tail file */
        sprintf(command,"cat %s/tail.txt | sed 's/latitu/%f/' > newTail.txt", 
                path, narr_lat[i]); 
        status = system(command);
        if (status != SUCCESS)
        {
            sprintf (errstr, "system call 2");
            LST_ERROR (errstr, "first_file");
        }

        sprintf(command,"cat newTail.txt | sed 's/longit/%f/' > newTail2.txt",
                    narr_lon[i]);
        status = system(command);
        if (status != SUCCESS)
        {
            sprintf (errstr, "system call 3");
            LST_ERROR (errstr, "first_file");
        }
 
        /* insert current julian day into tail file */
        sprintf(command,"cat newTail2.txt | sed 's/jay/%d/' > newTail3.txt",
                input->meta.acq_date.doy);
        status = system(command);
        if (status != SUCCESS)
        {
            sprintf (errstr, "system call 4");
            LST_ERROR (errstr, "first_file");
        }

        /* Allocate memory */
        temp_height = (float *)malloc(MAX_MODTRAN_LAYER * sizeof(float)); 
        if (temp_height == NULL)
        {
            sprintf (errstr, "Allocating temp_height memory");
            LST_ERROR (errstr, "first_files");
        }

        temp_pressure = (float *)malloc(MAX_MODTRAN_LAYER * sizeof(float)); 
        if (temp_pressure == NULL)
        {
            sprintf (errstr, "Allocating temp_pressure memory");
            LST_ERROR (errstr, "first_files");
        }

        temp_temp = (float *)malloc(MAX_MODTRAN_LAYER * sizeof(float)); 
        if (temp_temp == NULL)
        {
            sprintf (errstr, "Allocating temp_temp memory");
            LST_ERROR (errstr, "first_files");
        }

        temp_rh = (float *)malloc(MAX_MODTRAN_LAYER * sizeof(float)); 
        if (temp_rh == NULL)
        {
            sprintf (errstr, "Allocating temp_rh memory");
            LST_ERROR (errstr, "first_files");
        }

        /* iterature through all ground altitudes at which modtran is run */
        for (j = 0; j < NUM_ELEVATIONS; j++)
        {
            /* create a directory for the current height */
            sprintf(current_gdalt, "%s/%5.3f", current_point, gndalt[j]);

            if (stat(current_gdalt, &s)== -1)
            {
                status = mkdir(current_gdalt, 0755);
                if (status != SUCCESS)
                {
                    sprintf (errstr, "system call 4");
                    LST_ERROR (errstr, "first_file");
                }
            }

            /* determine layers below current gndalt and closest index 
               above and below */
            for (k = 0; k < NUM_ELEVATIONS; k++)
            {
                if ((narr_height[k][i] - gndalt[j]) >= MINSIGMA)
                {
                    index_below = k - 1;
                    index_above = k;
                    break;
                }
            }
                
            /* linearly interpolate pressure, temperature, and relative 
               humidity to gndalt for lowest layer */
            new_pressure = pressure[index_below][i] + (gndalt[j] -
                narr_height[index_below][i]) * ((pressure[index_above][i] - 
                pressure[index_below][i]) / (narr_height[index_above][i] - 
                narr_height[index_below][i]));
            new_temp = narr_tmp[index_below][i] + (gndalt[j] -
                narr_height[index_below][i]) * ((narr_tmp[index_above][i] - 
                narr_tmp[index_below][i]) / (narr_height[index_above][i] - 
                narr_height[index_below][i]));
            new_rh = narr_rh[index_below][i] + (gndalt[j] -
                narr_height[index_below][i]) * ((narr_rh[index_above][i] - 
                narr_rh[index_below][i]) / (narr_height[index_above][i] - 
                narr_height[index_below][i]));           

            /* create arrays containing only layers to be included in current 
               tape5 file */
            index = 0;
            if ((gndalt[j] - narr_height[index_above][i] - 0.001) < MINSIGMA)
            {
                for (k = index_above; k < NUM_ELEVATIONS; k++)
                {
                    temp_height[index] = narr_height[k][i];
                    temp_pressure[index] = pressure[k][i];
                    temp_temp[index] = narr_tmp[k][i];
                    temp_rh[index] = narr_rh[k][i];
                    index++;
                }
            }
            else
            {
                temp_height[index] = gndalt[j];
                temp_pressure[index] = new_pressure;
                temp_temp[index] = new_temp;
                temp_rh[index] = new_rh;
                index++;

                for (k = index_above; k < NUM_ELEVATIONS; k++)
                {
                    temp_height[index] = narr_height[k][i];
                    temp_pressure[index] = pressure[k][i];
                    temp_temp[index] = narr_tmp[k][i];
                    temp_rh[index] = narr_rh[k][i];
                    index++;
                }
            }

            /* determine maximum height of NARR layers and where the standard 
               atmosphere is greater than this */
            index2 = 0;
            for (k = 0; k < STAN_LAYER; k++)
            {
                if (stan_height[k] > narr_height[P_LAYER-1][i])
                {
                    counter[index2] = k;
                    index2++;
                }
            }

            /* if there are at least three layers above to highest NARR layer, 
               add standard atmosphere layers and linearly interpolate height, 
               pressure, temp, and rel hum to create a smooth transition between
               the NARR layers and the standard upper atmosphere */
            if (index2 >= 3)
            {
                new_height = (stan_height[counter[2]]+temp_height[index-1]) / 2.0; 
                new_pressure = temp_pressure[index-1]+(new_height-temp_height[index -1]) *
                    ((stan_pre[counter[2]]-temp_pressure[index - 1]) / 
                     (stan_height[counter[2]] - temp_height[index-1]));
                new_temp = temp_temp[index-1]+(new_height-temp_height[index-1]) * 
                    ((stan_temp[counter[2]]-temp_temp[index - 1]) / 
                    (stan_height[counter[2]] - temp_height[index-1]));            
                new_rh = temp_rh[counter[2]]+(new_height-temp_height[index-1]) * 
                    ((stan_rh[counter[2]]-temp_rh[index-1]) / 
                    (stan_height[counter[2]]-temp_height[index-1]));
            }


            /* concatenate NARR layers, new layer, and standard atmosphere layers */
            temp_height[index] = new_height;
            temp_pressure[index] = new_pressure;
            temp_temp[index] = new_temp;
            temp_rh[index] = new_rh;
            index++;
            for (k = 0; k < index2; k++)
            {
                temp_height[index] = stan_height[counter[k]];
                temp_pressure[index] = stan_pre[counter[k]];
                temp_temp[index] = stan_temp[counter[k]];
                temp_rh[index] = stan_rh[counter[k]];
                index++;
            }

            /* write atmospheric layers to a text file in format proper for tape5 file */
            fd = fopen("tempLayers.txt", "w"); 
            if (fd == NULL)
            {
                sprintf (errstr, "Opening file: tempLayers.txt\n");
                LST_ERROR (errstr, "first_files");
            }

            /* Write out the intermediate file */
            for (k = 0; k < index; k++)
            {
                fprintf(fd, "%10.3f,%10.3e,%10.3e,%10.3e,%10.3e,%10.3e,%16s\n", 
                        temp_height[k], temp_pressure[k], temp_temp[k],
                        temp_rh[k], 0.0, 0.0, "AAH             ");
            }

            /* Close the intermediate file */
            status = fclose(fd);
            if ( status )
            {
                sprintf (errstr, "Closing file: tempLayers.txt\n");
                LST_ERROR (errstr, "first_files");
            }

            /* determine number of layers for current ground altitude and insert 
               into head file */
            memcpy(temp_out, &index, sizeof(index)); 
            sprintf(command,"cat %s/head.txt | sed 's/nml/%s/' > newHead.txt", 
                    path, temp_out); 
            status = system(command);
            if (status != SUCCESS)
            {
                sprintf (errstr, "system call 5");
                LST_ERROR (errstr, "first_file");
            }

            /* insert current ground altitude into head file */
            memcpy(temp_out, &gndalt[j], sizeof(gndalt[j])); 
            sprintf(command,"cat newHead.txt | sed 's/gdalt/%s/' > newHead2.txt", temp_out); 
            status = system(command);
            if (status != SUCCESS)
            {
                sprintf (errstr, "system call 6");
                LST_ERROR (errstr, "first_file");
            }

            /* iterate through [temperature,albedo]  pairs at which to run modtran */
            for (k = 0; k <= 2; k++)
            {
                if (k == 2)
                    sprintf(temp_out, "%s", "000");
                else
                    sprintf(temp_out, "%d", (int) tmp[k]);

                /* create directory for the current temperature */
                sprintf(current_temp, "%s/%s", current_gdalt, temp_out);

                if (stat(current_temp, &s)== -1)
                {
                    status = mkdir(current_temp, 0755);
                    if (status != SUCCESS)
                    {
                        sprintf (errstr, "system call 7");
                        LST_ERROR (errstr, "first_file");
                    }
                }

                /* insert current temperature into head file */
                sprintf(command,"cat newHead2.txt | sed 's/tmp/%s/' > newHead3.txt", temp_out); 
                status = system(command);
                if (status != SUCCESS)
                {
                    sprintf (errstr, "system call 8");
                    LST_ERROR (errstr, "first_file");
                }

                /* create directory for the current albedo */
                sprintf(current_alb, "%s/%3.1f", current_temp, alb[k]);
                if (stat(current_alb, &s)== -1)
                {
                    status = mkdir(current_alb, 0755);
                    if (status != SUCCESS)
                    {
                        sprintf (errstr, "system call 9");
                        LST_ERROR (errstr, "first_file");
                    }
                }

                /* insert current albedo into head file */
                sprintf(command,"cat newHead3.txt | sed 's/alb/%f/' > newHead4.txt", alb[k]); 
                status = system(command);
                if (status != SUCCESS)
                {
                    sprintf (errstr, "system call 10");
                    LST_ERROR (errstr, "first_file");
                }
                
                /* concatenate head file, atmospheric layers, and tail file to create a 
                   tape5 file for modtran specific to this location and ground altitude 
                   with variables for temperature and albedo */
                sprintf(command, "cat newHead4.txt newTail3.txt tempLayers.txt > %s/tape5", 
                        current_alb); 
                status = system(command);
                if (status != SUCCESS)
                {
                    sprintf (errstr, "system call 11");
                    LST_ERROR (errstr, "first_file");
                }

                /* create string for case list containing location of current tape5 file
                   create string for command list containing commands for modtran run
                   iterate entry count*/
                *entry = 0;
                sprintf(case_list[*entry], "%s", current_alb);
                sprintf(command_list[*entry], "pushd %s ; ln -s %s/DATA; %s/Mod90_5.2.2.exe; popd'", 
                        modtran_path, modtran_path, current_alb); 
                entry++;
            }
        }
    }

    /* Free memory allocation */
    status = free_2d_array((void **)pressure);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: pressure\n");
        LST_ERROR (errstr, "first_files");              
    }

    status = free_2d_array((void **)narr_height);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: narr_height\n");
        LST_ERROR (errstr, "first_files");              
    }

    status = free_2d_array((void **)narr_rh);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: narr_rh\n");
        LST_ERROR (errstr, "first_files");              
    }

    status = free_2d_array((void **)narr_tmp);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: narr_tmp\n");
        LST_ERROR (errstr, "first_files");              
    }

    /* write caseList to a file */ 
    fd = fopen("caseList", "w"); 
    if (fd == NULL)
    {
        sprintf (errstr, "Opening file: caseList\n");
        LST_ERROR (errstr, "first_files");
    }

    /* Write out the caseList file */
    for (k = 0; k < *entry; k++)
    {
        fprintf(fd, "%s\n", case_list[k]);
    }

    /* Close the caseList file */
    status = fclose(fd);
    if ( status )
    {
        sprintf (errstr, "Closing file: caseList\n");
        LST_ERROR (errstr, "first_files");
    }

    /* write commandList to a file */
    fd = fopen("commandList", "w"); 
    if (fd == NULL)
    {
        sprintf (errstr, "Opening file: commandList\n");
        LST_ERROR (errstr, "first_files");
    }

    /* Write out the commandList file */
    for (k = 0; k < *entry; k++)
    {
        fprintf(fd, "%s\n", command_list[k]);
    }

    /* Close the commandList file */
    status = fclose(fd);
    if ( status )
    {
        sprintf (errstr, "Closing file: commandList\n");
        LST_ERROR (errstr, "first_files");
    }

    return SUCCESS;
}
