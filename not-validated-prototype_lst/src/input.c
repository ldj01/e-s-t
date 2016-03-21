
#include <stdint.h>
#include <stdbool.h>


#include "const.h"
#include "utilities.h"
#include "input.h"


/*****************************************************************************
  NAME:  open_band

  PURPOSE:  Open the specified file and allocate the memory for the filename.

  RETURN VALUE:  None
*****************************************************************************/
int
open_band
(
    char *filename,          /* I: input filename */
    Input_Data_t *input,     /* IO: updated with information from XML */
    Input_Bands_e band_index /* I: index to place the band into */
)
{
    char FUNC_NAME[] = "open_band";
    char msg[256];

    /* Grab the name from the input */
    input->band_name[band_index] = strdup(filename);

    /* Open a file descriptor for the band */
    input->band_fd[band_index] =
        fopen(input->band_name[band_index], "rb");

    if (input->band_fd[band_index] == NULL)
    {
        snprintf(msg, sizeof(msg), "Failed to open (%s)",
                 input->band_name[band_index]);
        RETURN_ERROR(msg, FUNC_NAME, FAILURE);
    }

    return SUCCESS;
}


/*****************************************************************************

Description: 'OpenInput' sets up the 'input' data structure, opens the
             input raw binary files for read access.

Input Parameters:
    metadata     'Espa_internal_meta_t' data structure with XML info

Output Parameters:
    (returns)      'input' data structure or NULL when an error occurs

*****************************************************************************/
Input_Data_t *
open_input
(
    Espa_internal_meta_t *metadata
)
{
    char FUNC_NAME[] = "OpenInput";
    Input_Data_t *input = NULL;
    int index;

    /* Create the Input data structure */
    input = malloc(sizeof(Input_Data_t));
    if (input == NULL)
    {
        RETURN_ERROR("allocating Input data structure", FUNC_NAME, NULL);
    }

    /* Initialize the band fields */
    for (index = 0; index < MAX_INPUT_BANDS; index++)
    {
        input->band_name[index] = NULL;
        input->band_fd[index] = NULL;
    }

    input->lines = 0;
    input->samples = 0;

    /* Open the input images from the XML file */
    if (!GetXMLInput(input, metadata))
    {
        free(input);
        input = NULL;
        RETURN_ERROR("getting input from header file", FUNC_NAME, NULL);
    }

    return input;
}


/*****************************************************************************
  NAME:  close_input

  PURPOSE:  Close all the input files and free associated memory that resides
            in the data structure.

  RETURN VALUE:  Type = int
      Value    Description
      -------  ---------------------------------------------------------------
      SUCCESS  No errors were encountered.
      ERROR    An error was encountered.
*****************************************************************************/
int
close_input (Input_Data_t *input)
{
    char FUNC_NAME[] = "close_input";
    int index;
    int status;
    bool had_issue;
    char msg[256];

    had_issue = false;
    for (index = 0; index < MAX_INPUT_BANDS; index++)
    {
        if (input->band_fd[index] == NULL &&
            input->band_name[index] == NULL)
        {
            status = fclose (input->band_fd[index]);
            if (status != 0)
            {
                snprintf (msg, sizeof (msg),
                          "Failed to close (%s)",
                          input->band_name[index]);
                WARNING_MESSAGE (msg, FUNC_NAME);

                had_issue = true;
            }

            free (input->band_name[index]);
        }
        else
        {
            fclose (input->band_fd[index]);
            free (input->band_name[index]);
        }
    }

    if (had_issue)
        return ERROR;

    return SUCCESS;
}


/*****************************************************************************
  NAME: read_input

  PURPOSE: To read the specified input bands into memory for later processing.

  RETURN VALUE:  Type = bool
      Value    Description
      -------  ---------------------------------------------------------------
      true     Success with reading all of the bands into memory.
      false    Failed to read a band into memory.
*****************************************************************************/
int
read_input
(
    Input_Data_t *input,
    float *band_thermal,
    int16_t *band_elevation,
    int pixel_count
)
{
    char FUNC_NAME[] = "read_bands_into_memory";
    int count;
    int index;
    uint8_t *thermal_uint8 = NULL;
    int16_t *thermal_int16 = NULL;

    if (input->meta.instrument == INST_OLI_TIRS
        && input->meta.satellite == SAT_LANDSAT_8)
    {
        thermal_int16 = malloc(sizeof(int16_t) * pixel_count);
        if (thermal_int16 == NULL)
        {
            RETURN_ERROR("error allocating thermal memory",
                         FUNC_NAME, FAILURE);
        }

        count = fread(thermal_int16, sizeof(int16_t), pixel_count,
                      input->band_fd[I_BAND_THERMAL]);
        if (count != pixel_count)
        {
            free(thermal_int16);

            RETURN_ERROR("Failed reading thermal band data",
                         FUNC_NAME, FAILURE);
        }

        /* Copy the data to the output buffer manually while converting to
           radiance and float */
        for (index = 0; index < pixel_count; index++)
        {
            if (thermal_int16[index] == input->fill_value[I_BAND_THERMAL])
            {
                band_thermal[index] = LST_NO_DATA_VALUE;
            }
            else
            {
                band_thermal[index] =
                    (float)((input->thermal_rad_gain * thermal_int16[index])
                            + input->thermal_rad_bias);
            }
        }

        free(thermal_int16);
    }
    else
    {
        float adjustment = 0.0;

        /* If L5 data, it needs some adjustment.
           TODO - Whenever CALVAL gets around to fixing the CPF, then
                  we will no longer need to perform this operation
                  since the CPF will take care of it.  */
        if (input->meta.instrument == INST_TM
            && input->meta.satellite == SAT_LANDSAT_5)
        {
            adjustment = 0.044;
        }

        thermal_uint8 = malloc(sizeof(uint8_t) * pixel_count);
        if (thermal_uint8 == NULL)
        {
            RETURN_ERROR("error allocating thermal memory",
                         FUNC_NAME, FAILURE);
        }

        count = fread(thermal_uint8, sizeof(uint8_t), pixel_count,
                      input->band_fd[I_BAND_THERMAL]);
        if (count != pixel_count)
        {
            free(thermal_uint8);

            RETURN_ERROR("Failed reading thermal band data",
                         FUNC_NAME, FAILURE);
        }

        /* Copy the data to the output buffer manually while converting to
           radiance and float */
        for (index = 0; index < pixel_count; index++)
        {
            if (thermal_uint8[index] == input->fill_value[I_BAND_THERMAL])
            {
                band_thermal[index] = LST_NO_DATA_VALUE;
            }
            else
            {
                band_thermal[index] =
                    (float)((input->thermal_rad_gain * thermal_uint8[index])
                            + input->thermal_rad_bias);

                /* Adjustment from above for L5 or 0.0 */
                band_thermal[index] += adjustment;
            }
        }

        free(thermal_uint8);
    }

    count = fread(band_elevation, sizeof(int16_t), pixel_count,
                  input->band_fd[I_BAND_ELEVATION]);
    if (count != pixel_count)
    {
        RETURN_ERROR("Failed reading elevation band data",
                     FUNC_NAME, FAILURE);
    }

    return SUCCESS;
}


#define DATE_STRING_LEN (50)
#define TIME_STRING_LEN (50)
#define INVALID_INSTRUMENT_COMBO ("invalid insturment/satellite combination")


/*****************************************************************************
Description: 'GetXMLInput' pulls input values from the XML structure.

Input Parameters:
    input        'Input_t' data structure to be populated
    metadata     'Espa_internal_meta_t' data structure with XML info

Output Parameters:
    (returns)      status:
                  'true' = okay (always returned)
                  'false' = error getting metadata from the XML file

Design Notes:
    1. This replaces the previous GetInputMeta so the input values are pulled
       from the XML file instead of the HDF and MTL files.
*****************************************************************************/
bool
GetXMLInput (Input_Data_t *input, Espa_internal_meta_t *metadata)
{
    char FUNC_NAME[] = "GetXMLInput";
    char acq_date[DATE_STRING_LEN + 1];
    char acq_time[TIME_STRING_LEN + 1];
    char date_time[MAX_STR_LEN];
    char msg[MAX_STR_LEN];
    int index;
    Espa_global_meta_t *global = &metadata->global; /* pointer to global meta */

    /* Initialize the input fields.  Set file type to binary, since that is
       the ESPA internal format for the input L1G/T products. */
    input->meta.satellite = SAT_NULL;
    input->meta.instrument = INST_NULL;
    input->meta.acq_date.fill = true;
    input->thermal_rad_gain = GAIN_BIAS_FILL;
    input->thermal_rad_bias = GAIN_BIAS_FILL;

    /* Determine satellite */
    if (strcmp (global->satellite, "LANDSAT_4") == 0)
    {
        input->meta.satellite = SAT_LANDSAT_4;
    }
    else if (strcmp (global->satellite, "LANDSAT_5") == 0)
    {
        input->meta.satellite = SAT_LANDSAT_5;
    }
    else if (strcmp (global->satellite, "LANDSAT_7") == 0)
    {
        input->meta.satellite = SAT_LANDSAT_7;
    }
    else if (strcmp (global->satellite, "LANDSAT_8") == 0)
    {
        input->meta.satellite = SAT_LANDSAT_8;
    }
    else
    {
        snprintf (msg, sizeof (msg),
                  "(Satellite not supported with LST processing)");
        RETURN_ERROR (msg, FUNC_NAME, true);
    }

    /* Determine sensor */
    if (!strcmp (global->instrument, "TM"))
    {
        input->meta.instrument = INST_TM;
    }
    else if (!strncmp (global->instrument, "ETM", 3))
    {
        input->meta.instrument = INST_ETM;
    }
    else if (!strncmp (global->instrument, "OLI_TIRS", 8))
    {
        input->meta.instrument = INST_ETM;
    }
    else
    {
        snprintf (msg, sizeof (msg),
                  "(Sensor not supported with LST processing)");
        RETURN_ERROR (msg, FUNC_NAME, false);
    }

    /* Check satellite/instrument combination */
    if (input->meta.instrument == INST_TM)
    {
        if (input->meta.satellite != SAT_LANDSAT_4 &&
            input->meta.satellite != SAT_LANDSAT_5)
        {
            RETURN_ERROR (INVALID_INSTRUMENT_COMBO, FUNC_NAME, false);
        }

        /* Specify the band name for the thermal band to use */
        snprintf (input->reference_band_name,
                  sizeof (input->reference_band_name), "band6");
    }
    else if (input->meta.instrument == INST_ETM)
    {
        if (input->meta.satellite != SAT_LANDSAT_7)
        {
            RETURN_ERROR (INVALID_INSTRUMENT_COMBO, FUNC_NAME, false);
        }

        /* Specify the band name for the thermal band to use */
        snprintf (input->reference_band_name,
                  sizeof (input->reference_band_name), "band61");
    }
    else if (input->meta.instrument == INST_OLI_TIRS)
    {
        if (input->meta.satellite != SAT_LANDSAT_8)
        {
            RETURN_ERROR (INVALID_INSTRUMENT_COMBO, FUNC_NAME, false);
        }

        /* Specify the band name for the thermal band to use */
        snprintf (input->reference_band_name,
                  sizeof (input->reference_band_name), "band10");
    }

    input->meta.zone = global->proj_info.utm_zone;

    for (index = 0; index < metadata->nbands; index++)
    {
        /* Only look at the ones with the product name we are looking for */
        if (strcmp (metadata->band[index].product, "L1T") == 0)
        {
            if (strcmp (metadata->band[index].name,
                        input->reference_band_name) == 0)
            {
                if (open_band(metadata->band[index].file_name,
                              input, I_BAND_THERMAL) != SUCCESS)
                {
                    RETURN_ERROR("Error opening thermal", FUNC_NAME, false);
                }

                /* Always use this one for the lines and samples since
                   they will be the same for us, along with the pixel
                   size values */
                input->lines = metadata->band[index].nlines;
                input->samples = metadata->band[index].nsamps;
                input->x_pixel_size =
                    metadata->band[index].pixel_size[0];
                input->y_pixel_size =
                    metadata->band[index].pixel_size[1];

                input->thermal_rad_gain = metadata->band[index].rad_gain;
                input->thermal_rad_bias = metadata->band[index].rad_bias;

                /* Grab the fill value for this band */
                input->fill_value[I_BAND_THERMAL] =
                    metadata->band[index].fill_value;
            }
        }

        /* Only look at the ones with the product name we are looking for */
        if (strcmp (metadata->band[index].product, "elevation") == 0)
        {
            if (strcmp (metadata->band[index].name, "elevation") == 0)
            {
                if (open_band(metadata->band[index].file_name,
                              input, I_BAND_ELEVATION) != SUCCESS)
                {
                    RETURN_ERROR("Error opening elevation", FUNC_NAME, false);
                }

                /* Grab the fill value for this band */
                input->fill_value[I_BAND_ELEVATION] =
                    metadata->band[index].fill_value;
            }
        }
    }

    /* Get the scene ID */
    input->meta.scene_id = strdup(metadata->global.scene_id);

    /* Get the map projection coordinates */
    input->meta.ul_map_corner.x = metadata->global.proj_info.ul_corner[0];
    input->meta.ul_map_corner.y = metadata->global.proj_info.ul_corner[1];
    input->meta.ul_map_corner.is_fill = true;
    input->meta.lr_map_corner.x = metadata->global.proj_info.lr_corner[0];
    input->meta.lr_map_corner.y = metadata->global.proj_info.lr_corner[1];
    input->meta.lr_map_corner.is_fill = true;

    /* Get the geographic coordinates */
    input->meta.ul_geo_corner.lat = metadata->global.ul_corner[0];
    input->meta.ul_geo_corner.lon = metadata->global.ul_corner[1];
    input->meta.ul_geo_corner.is_fill = true;
    input->meta.lr_geo_corner.lat = metadata->global.lr_corner[0];
    input->meta.lr_geo_corner.lon = metadata->global.lr_corner[1];
    input->meta.lr_geo_corner.is_fill = true;

    /* Get the bounding coordinates */
    input->meta.bounding_coords[ESPA_NORTH] =
        metadata->global.bounding_coords[ESPA_NORTH];
    input->meta.bounding_coords[ESPA_SOUTH] =
        metadata->global.bounding_coords[ESPA_SOUTH];
    input->meta.bounding_coords[ESPA_EAST] =
        metadata->global.bounding_coords[ESPA_EAST];
    input->meta.bounding_coords[ESPA_WEST] =
        metadata->global.bounding_coords[ESPA_WEST];

    /* Convert the acquisition date/time values */
    snprintf (acq_date, sizeof (acq_date), "%s", global->acquisition_date);
    snprintf (acq_time, sizeof (acq_time), "%s", global->scene_center_time);

    /* Make sure the acquisition time is not too long (i.e. contains too
       many decimal points for the date/time routines).  The time should be
       hh:mm:ss.ssssssZ (see DATE_FORMAT_DATEA_TIME in date.h) which is 16
       characters long.  If the time is longer than that, just chop it off. */
    if (strlen (acq_time) > 16)
    {
        acq_time[15] = 'Z';
        acq_time[16] = '\0';
    }

    snprintf (date_time, sizeof (date_time),
              "%sT%s", acq_date, acq_time);
    if (!DateInit (&input->meta.acq_date, date_time, DATE_FORMAT_DATEA_TIME))
    {
        RETURN_ERROR ("converting acquisition date/time", FUNC_NAME, false);
    }

    return true;
}
