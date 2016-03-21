
#include "const.h"
#include "input.h"
#include "utilities.h"
#include "intermediate_data.h"


int
open_intermediate(Input_Data_t *input,
                  Intermediate_Data_t *inter)
{
    char *FUNC_NAME = "open_intermediate";
    char msg[PATH_MAX];

    /* First figure out and assign the filenames */
    snprintf(inter->thermal_filename,
             sizeof(inter->thermal_filename),
             "%s_%s.img",
             input->meta.scene_id,
             LST_THERMAL_RADIANCE_BAND_NAME);
    snprintf(inter->upwelled_filename,
             sizeof(inter->upwelled_filename),
             "%s_%s.img",
             input->meta.scene_id,
             LST_UPWELLED_RADIANCE_BAND_NAME);
    snprintf(inter->downwelled_filename,
             sizeof(inter->downwelled_filename),
             "%s_%s.img",
             input->meta.scene_id,
             LST_DOWNWELLED_RADIANCE_BAND_NAME);
    snprintf(inter->transmittance_filename,
             sizeof(inter->transmittance_filename),
             "%s_%s.img",
             input->meta.scene_id,
             LST_ATMOS_TRANS_BAND_NAME);

    /* Now open the file descriptors */
    inter->thermal_fd = fopen(inter->thermal_filename, "wb");
    if (inter->thermal_fd == NULL)
    {
        sprintf(msg, "Opening intermediate file: %s",
                inter->thermal_filename);
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }

    inter->transmittance_fd = fopen(inter->transmittance_filename, "wb");
    if (inter->transmittance_fd == NULL)
    {
        sprintf(msg, "Opening intermediate file: %s",
                inter->transmittance_filename);
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }

    inter->upwelled_fd = fopen(inter->upwelled_filename, "wb");
    if (inter->upwelled_fd == NULL)
    {
        sprintf(msg, "Opening intermediate file: %s",
                inter->upwelled_filename);
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }

    inter->downwelled_fd = fopen(inter->downwelled_filename, "wb");
    if (inter->downwelled_fd == NULL)
    {
        sprintf(msg, "Opening intermediate file: %s",
                inter->downwelled_filename);
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }

    /* Initialize the memory items */
    inter->band_thermal = NULL;
    inter->band_transmittance = NULL;
    inter->band_upwelled = NULL;
    inter->band_downwelled = NULL;

#if OUTPUT_CELL_DESIGNATION_BAND
    snprintf(inter->cell_filename,
             sizeof(inter->cell_filename),
             "%s_cellnumbers.img",
             input->meta.scene_id);

    inter->cell_fd = fopen(inter->cell_filename, "wb");
    if (inter->cell_fd == NULL)
    {
        sprintf(msg, "Opening intermediate file: %s",
                inter->cell_filename);
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }

    inter->band_cell = NULL;
#endif

    return SUCCESS;
}


int
write_intermediate(Intermediate_Data_t *inter,
                   int pixel_count)
{
    char *FUNC_NAME = "close_intermediate";
    char msg[PATH_MAX];
    int status;

    status = fwrite(inter->band_thermal,
                    sizeof(float),
                    pixel_count,
                    inter->thermal_fd);
    if (status != pixel_count)
    {
        sprintf (msg, "Writing to %s", inter->thermal_filename);
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }

    status = fwrite(inter->band_transmittance,
                    sizeof(float),
                    pixel_count,
                    inter->transmittance_fd);
    if (status != pixel_count)
    {
        sprintf (msg, "Writing to %s", inter->transmittance_filename);
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }

    status = fwrite(inter->band_upwelled,
                    sizeof(float),
                    pixel_count,
                    inter->upwelled_fd);
    if (status != pixel_count)
    {
        sprintf (msg, "Writing to %s", inter->upwelled_filename);
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }

    status = fwrite(inter->band_downwelled,
                    sizeof(float),
                    pixel_count,
                    inter->downwelled_fd);
    if (status != pixel_count)
    {
        sprintf (msg, "Writing to %s", inter->downwelled_filename);
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }

#if OUTPUT_CELL_DESIGNATION_BAND
    status = fwrite(inter->band_cell,
                    sizeof(uint8_t),
                    pixel_count,
                    inter->cell_fd);
    if (status != pixel_count)
    {
        sprintf (msg, "Writing to %s", inter->cell_filename);
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }
#endif

    return SUCCESS;
}


int
close_intermediate(Intermediate_Data_t *inter)
{
    char *FUNC_NAME = "close_intermediate";
    char msg[PATH_MAX];
    int status;

    status = fclose(inter->thermal_fd);
    if (status)
    {
        sprintf(msg, "Closing file %s", inter->thermal_filename);
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }
    inter->thermal_fd = NULL;

    status = fclose(inter->transmittance_fd);
    if (status)
    {
        sprintf(msg, "Closing file %s", inter->transmittance_filename);
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }
    inter->transmittance_fd = NULL;

    status = fclose(inter->upwelled_fd);
    if (status)
    {
        sprintf(msg, "Closing file %s", inter->upwelled_filename);
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }
    inter->upwelled_fd = NULL;

    status = fclose(inter->downwelled_fd);
    if (status)
    {
        sprintf(msg, "Closing file %s", inter->downwelled_filename);
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }
    inter->downwelled_fd = NULL;

#if OUTPUT_CELL_DESIGNATION_BAND
    status = fclose(inter->cell_fd);
    if (status)
    {
        sprintf(msg, "Closing file %s", inter->cell_filename);
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }
    inter->cell_fd = NULL;
#endif

    return SUCCESS;
}


int
allocate_intermediate(Intermediate_Data_t *inter,
                      int pixel_count)
{
    char *FUNC_NAME = "allocate_intermediate";
    char msg[PATH_MAX];

    inter->band_thermal = calloc(pixel_count, sizeof(float));
    if (inter->band_thermal == NULL)
    {
        sprintf(msg, "Allocating memory for %s",
                inter->thermal_filename);
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }

    inter->band_transmittance = calloc(pixel_count, sizeof(float));
    if (inter->band_transmittance == NULL)
    {
        free_intermediate(inter);

        sprintf(msg, "Allocating memory for %s",
                inter->transmittance_filename);
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }

    inter->band_upwelled = calloc(pixel_count, sizeof(float));
    if (inter->band_upwelled == NULL)
    {
        free_intermediate(inter);

        sprintf(msg, "Allocating memory for %s",
                inter->upwelled_filename);
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }

    inter->band_downwelled = calloc(pixel_count, sizeof(float));
    if (inter->band_downwelled == NULL)
    {
        free_intermediate(inter);

        sprintf(msg, "Allocating memory for %s",
                inter->downwelled_filename);
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }

#if OUTPUT_CELL_DESIGNATION_BAND
    inter->band_cell = calloc(pixel_count, sizeof(uint8_t));
    if (inter->band_cell == NULL)
    {
        free_intermediate(inter);

        sprintf(msg, "Allocating memory for %s",
                inter->cell_filename);
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }
#endif

    return SUCCESS;
}


void
free_intermediate(Intermediate_Data_t *inter)
{
    free(inter->band_thermal);
    inter->band_thermal = NULL;

    free(inter->band_transmittance);
    inter->band_transmittance = NULL;

    free(inter->band_upwelled);
    inter->band_upwelled = NULL;

    free(inter->band_downwelled);
    inter->band_downwelled = NULL;

#if OUTPUT_CELL_DESIGNATION_BAND
    free(inter->band_cell);
    inter->band_cell = NULL;
#endif
}
