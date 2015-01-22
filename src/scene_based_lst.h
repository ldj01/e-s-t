#ifndef SCENE_BASED_LST_H
#define SCENE_BASED_LST_H


#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>


#include "input.h"

typedef struct
{
    char full_path[PATH_MAX];
    float latitude;
    float longitude;
    float height;
} CASE_POINT;


int get_args
(
    int argc,                /* I: number of cmd-line args */
    char *argv[],            /* I: string of cmd-line args */
    char *xml_infile,        /* I: address of input XML metadata filename  */
    char *dem_infile,        /* I: address of input DEM filename */
    char *emissivity_infile, /* I: address of input emissivity filename */
    bool *tape_6,            /* O: use the tape6 output */
    bool *verbose,           /* O: verbose flag */
    bool *debug              /* O: debug flag */
);


int build_modtran_input
(
    Input_t *input,       /* I: input structure */
    int *num_pts,         /* O: number of NARR points */
    int *num_runs,        /* O: number of MODTRAN runs */
    CASE_POINT **case_list, /* O: case list information (allocated here) */
    char ***command_list, /* O: command list information (allocated here) */
    bool verbose,         /* I: value to indicate if intermediate messages
                                will be printed */
    bool debug            /* I: value to indicate if debug should be
                                generated */
);


int second_narr
(
    Input_t *input,   /* I: input structure */
    int num_points,   /* I: number of narr points */
    float alb,        /* I: albedo */
    CASE_POINT *case_list, /* I: modtran run list */
    float **results,  /* O: atmospheric parameter for modtarn run */
    bool verbose      /* I: value to indicate if intermediate messages will be
                            printed */
);


int third_pixels_post
(
    Input_t *input,          /* I: input structure */
    int num_points,          /* I: number of narr points */
    char *dem_infile,        /* I: address of input DEM filename */
    char *emi_infile,        /* I: address of input Emissivity filename */
    float **modtran_results, /* I: atmospheric parameter for modtarn run */
    bool verbose             /* I: value to indicate if intermediate messages
                                   will be printed */
);


#endif /* SCENE_BASED_LST_H */
