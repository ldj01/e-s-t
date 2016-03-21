
#ifndef INTERMEDIATE_DATA_H
#define INTERMEDIATE_DATA_H


#include "input.h"

/* This is for compile time debugging logic.
   Set it to 0 to turn it off.
   It needs to be set to 0 for production/standard processing. */
#define OUTPUT_CELL_DESIGNATION_BAND 1


/* Structure for the intermediate data */
typedef struct
{
    char thermal_filename[PATH_MAX];
    char transmittance_filename[PATH_MAX];
    char upwelled_filename[PATH_MAX];
    char downwelled_filename[PATH_MAX];
    FILE *thermal_fd;
    FILE *transmittance_fd;
    FILE *upwelled_fd;
    FILE *downwelled_fd;
    float *band_thermal;
    float *band_transmittance;
    float *band_upwelled;
    float *band_downwelled;

#if OUTPUT_CELL_DESIGNATION_BAND
    char cell_filename[PATH_MAX];
    FILE *cell_fd;
    uint8_t *band_cell;
#endif
} Intermediate_Data_t;


int open_intermediate(Input_Data_t *input,
                      Intermediate_Data_t *inter);

int write_intermediate(Intermediate_Data_t *inter,
                       int pixel_count);

int close_intermediate(Intermediate_Data_t *inter);

int allocate_intermediate(Intermediate_Data_t *inter,
                          int pixel_count);

void free_intermediate(Intermediate_Data_t *inter);


#endif /* INTERMEDIATE_DATA_H */
