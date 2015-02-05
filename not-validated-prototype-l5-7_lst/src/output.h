#ifndef OUTPUT_H
#define OUTPUT_H


#include <time.h>


#include "input.h"


#define FILL_VALUE 255


/* Structure for the 'output' data type */

typedef struct
{
    bool open;             /* Flag to indicate whether output file is open 
                              for access; 'true' = open, 'false' = not open */
    Img_coord_int_t size;  /* Output image size */
    int nband;             /* Number of output image bands */
    Espa_internal_meta_t metadata; /* metadata container to hold the band
                                      metadata for the output band; global
                                      metadata won't be valid */
    FILE *fp_bin;          /* File pointer for binary output file */
} Output_t;


/* Prototypes */
Output_t *OpenOutput (Espa_internal_meta_t *in_meta, Input_t *input);


bool PutOutput (Output_t *this, unsigned char **final_mask);


bool CloseOutput (Output_t *this);


bool FreeOutput (Output_t *this);


#endif
