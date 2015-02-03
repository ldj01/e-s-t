#ifndef BUILD_MODTRAN_INPUT_H
#define BUILD_MODTRAN_INPUT_H


#include <stdbool.h>


#include "lst_types.h"
#include "input.h"


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


#endif /* BUILD_MODTRAN_INPUT_H */
