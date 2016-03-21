
#ifndef BUILD_MODTRAN_INPUT_H
#define BUILD_MODTRAN_INPUT_H


#include <stdbool.h>


#include "lst_types.h"
#include "input.h"


int build_modtran_input
(
    Input_Data_t *input,       /* I: input structure */
    REANALYSIS_POINTS *points, /* I/O: The coordinate points */
    bool verbose,         /* I: value to indicate if intermediate messages
                                will be printed */
    bool debug            /* I: value to indicate if debug should be
                                generated */
);


#endif /* BUILD_MODTRAN_INPUT_H */
