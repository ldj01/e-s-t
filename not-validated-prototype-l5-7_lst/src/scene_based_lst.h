#ifndef CFMASK_H
#define CFMASK_H
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>
#include <input.h>
#include "error.h"

int first_files
(
    Input_t *input,             /*I: input structure */
    char **case_list,           /*O: modtran run list */
    char **command_list,        /*O: modtran run command list */
    int *num_points,            /*O: number of NARR points */
    bool verbose                /*I: value to indicate if intermediate messages 
                                     be printed */
);

int second_narr
(
    Input_t *input,             /*I: input structure */
    int num_points,             /*I: number of narr points */
    float alb,                  /*I: albedo */ 
    char **case_list,           /*I: modtran run list */
    float **results,            /*O: atmospheric parameter for modtarn run */
    bool verbose                /*I: value to indicate if intermediate messages 
                                     be printed */
);

int third_pixels_post
(
    Input_t *input,             /*I: input structure */
    int num_points,             /*I: number of narr points */
    char *dem_infile,           /*I: address of input DEM filename */
    char *emi_infile,           /*I: address of input Emissivity filename */
    float **modtran_results,    /*I: atmospheric parameter for modtarn run */
    bool verbose                /*I: value to indicate if intermediate messages 
                                     be printed */
);

#endif
