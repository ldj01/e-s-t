
#include <stdint.h>
#include <stdbool.h>

#include "const.h"
#include "utilities.h"
#include "input.h"


/* Functions */
Input_t *
OpenInput
(
    Espa_internal_meta_t *metadata
)
/*****************************************************************************

Description: 'OpenInput' sets up the 'input' data structure, opens the
             input raw binary files for read access.

Input Parameters:
    metadata     'Espa_internal_meta_t' data structure with XML info

Output Parameters:
    (returns)      'input' data structure or NULL when an error occurs

Team Unique Header:

*****************************************************************************/
{
    char FUNC_NAME[] = "OpenInput";
    Input_t *input = NULL;

    /* Create the Input data structure */
    input = malloc (sizeof (Input_t));
    if (input == NULL)
    {
        RETURN_ERROR ("allocating Input data structure", FUNC_NAME, NULL);
    }

    /* Initialize and get input from header file */
    if (!GetXMLInput (input, metadata))
    {
        free (input);
        input = NULL;
        RETURN_ERROR ("getting input from header file", FUNC_NAME, NULL);
    }

    /* Open file for access */
    input->thermal.fd = fopen (input->thermal.filename, "r");
    if (input->thermal.fd == NULL)
    {
        free (input->thermal.filename);
        input->thermal.filename = NULL;

        fclose (input->thermal.fd);

        free (input);
        input = NULL;

        RETURN_ERROR ("opening thermal file", FUNC_NAME, NULL);
    }
    else
        input->thermal.is_open = true;

    return input;
}


bool
GetInputThermLine
(
    Input_t *input,
    int iline,
    float *thermal_data
)
{
    char FUNC_NAME[] = "GetInputThermLine";
    void *buf = NULL;
    long loc;         /* pointer location in the raw binary file */
    int sample;
    uint8_t *line_uint8 = NULL;
    int16_t *line_int16 = NULL;

    /* Check the parameters */
    if (input == NULL)
    {
        RETURN_ERROR ("invalid input structure", FUNC_NAME, false);
    }

    if (!input->thermal.is_open)
    {
        RETURN_ERROR ("file not open", FUNC_NAME, false);
    }

    if (iline < 0 || iline >= input->thermal.size.l)
    {
        RETURN_ERROR ("invalid line number", FUNC_NAME, false);
    }

    /* Read the data */
    if (input->meta.instrument == INST_OLI_TIRS
        && input->meta.satellite == SAT_LANDSAT_8)
    {
        line_int16 = malloc (input->thermal.size.s * sizeof (int16_t));
        if (line_int16 == NULL)
        {
            RETURN_ERROR ("error allocating memory", FUNC_NAME, false);
        }

        buf = (void *) thermal_data;
        loc = (long) (iline * input->thermal.size.s * sizeof (int16_t));
        if (fseek (input->thermal.fd, loc, SEEK_SET))
        {
            RETURN_ERROR ("error seeking thermal line (binary)",
                          FUNC_NAME, false);
        }

        if (read_raw_binary (input->thermal.fd, 1, input->thermal.size.s,
                             sizeof (int16_t), buf)
            != SUCCESS)
        {
            RETURN_ERROR ("error reading thermal line (binary)",
                          FUNC_NAME, false);
        }

        /* Copy the data to the output buffer manually while converting to
           radiance and float */
        for (sample = 0; sample < input->thermal.size.s; sample++)
        {
            if (line_int16[sample] == input->thermal.fill_value)
            {
                thermal_data[sample] = LST_FILL_VALUE;
            }
            else
            {
                thermal_data[sample] =
                    (float) ((input->thermal.toa_gain * line_int16[sample])
                             + input->thermal.toa_bias);
            }
        }

        free (line_int16);
    }
    else
    {
        line_uint8 = malloc (input->thermal.size.s * sizeof (uint8_t));
        if (line_uint8 == NULL)
        {
            RETURN_ERROR ("error allocating memory", FUNC_NAME, false);
        }

        buf = (void *) line_uint8;
        loc = (long) (iline * input->thermal.size.s * sizeof (uint8_t));
        if (fseek (input->thermal.fd, loc, SEEK_SET))
        {
            RETURN_ERROR ("error seeking thermal line (binary)",
                          FUNC_NAME, false);
        }

        if (read_raw_binary (input->thermal.fd, 1, input->thermal.size.s,
                             sizeof (uint8_t), buf)
            != SUCCESS)
        {
            RETURN_ERROR ("error reading thermal line (binary)",
                          FUNC_NAME, false);
        }

        /* Copy the data to the output buffer manually while converting to
           radiance and float */
        for (sample = 0; sample < input->thermal.size.s; sample++)
        {
            if (line_uint8[sample] == input->thermal.fill_value)
            {
                thermal_data[sample] = LST_FILL_VALUE;
            }
            else
            {
                thermal_data[sample] =
                    (float) ((input->thermal.toa_gain * line_uint8[sample])
                             + input->thermal.toa_bias);

                /* If L5 data, it needs some adjustment, I don't know why...
                   ???? it was in the original code ???? */
                if (input->meta.instrument == INST_TM
                    && input->meta.satellite == SAT_LANDSAT_5)
                {
                    thermal_data[sample] += 0.044;
                }
            }
        }

        free (line_uint8);
    }

    return true;
}


bool
CloseInput (Input_t *input)
/*****************************************************************************

Description: 'CloseInput' ends SDS access and closes the input file.

Input Parameters:
     input     'input' data structure

Output Parameters:
     input     'input' data structure; the following fields are modified:
               open
     (returns) status:
               'true' = okay
               'false' = error return

*****************************************************************************/
{
    char FUNC_NAME[] = "CloseInput";

    if (input == NULL)
    {
        RETURN_ERROR ("invalid input structure", FUNC_NAME, false);
    }

    /*** now close the thermal file ***/
    if (input->thermal.is_open)
    {
        fclose (input->thermal.fd);
        input->thermal.is_open = false;
    }

    return true;
}


bool
FreeInput (Input_t *input)
/*****************************************************************************

Description: 'FreeInput' frees the 'input' data structure memory.

Input Parameters:
    input          'input' data structure

Output Parameters:
    (returns)      status:
                  'true' = okay (always returned)

*****************************************************************************/
{
    free (input->thermal.filename);
    input->thermal.filename = NULL;

    free (input);
    input = NULL;

    return true;
}


#define DATE_STRING_LEN (50)
#define TIME_STRING_LEN (50)
#define INVALID_INSTRUMENT_COMBO ("invalid insturment/satellite combination")


bool
GetXMLInput (Input_t *input, Espa_internal_meta_t * metadata)
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
{
    char FUNC_NAME[] = "GetXMLInput";
    char acq_date[DATE_STRING_LEN + 1];
    char acq_time[TIME_STRING_LEN + 1];
    char date_time[MAX_STR_LEN];
    char band_name[30];
    char msg[MAX_STR_LEN];
    int index;
    Espa_global_meta_t *global = &metadata->global; /* pointer to global meta */

    /* Initialize the input fields.  Set file type to binary, since that is
       the ESPA internal format for the input L1G/T products. */
    input->meta.satellite = SAT_NULL;
    input->meta.instrument = INST_NULL;
    input->meta.acq_date.fill = true;
    input->thermal.filename = NULL;
    input->thermal.fd = NULL;
    input->thermal.is_open = false;
    input->thermal.size.s = -1;
    input->thermal.size.l = -1;
    input->thermal.toa_gain = GAIN_BIAS_FILL;
    input->thermal.toa_bias = GAIN_BIAS_FILL;
    input->thermal.fill_value = INPUT_FILL;

    /* Determine satellite */
    if (strcmp (global->satellite, "LANDSAT_5") == 0)
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
        RETURN_ERROR (msg, FUNC_NAME, true);
    }

    /* Check satellite/instrument combination */
    if (input->meta.instrument == INST_TM)
    {
        if (input->meta.satellite != SAT_LANDSAT_5)
        {
            RETURN_ERROR (INVALID_INSTRUMENT_COMBO, FUNC_NAME, true);
        }

        /* Specify the band name for the thermal band to use */
        snprintf (band_name, sizeof (band_name), "band6");
    }
    else if (input->meta.instrument == INST_ETM)
    {
        if (input->meta.satellite != SAT_LANDSAT_7)
        {
            RETURN_ERROR (INVALID_INSTRUMENT_COMBO, FUNC_NAME, true);
        }

        /* Specify the band name for the thermal band to use */
        snprintf (band_name, sizeof (band_name), "band61");
    }
    else if (input->meta.instrument == INST_OLI_TIRS)
    {
        if (input->meta.satellite != SAT_LANDSAT_8)
        {
            RETURN_ERROR (INVALID_INSTRUMENT_COMBO, FUNC_NAME, true);
        }

        /* Specify the band name for the thermal band to use */
        snprintf (band_name, sizeof (band_name), "band10");
    }

    input->meta.zone = global->proj_info.utm_zone;

    for (index = 0; index < metadata->nbands; index++)
    {
        /* Only look at the ones with the product name we are looking for */
        if (strcmp (metadata->band[index].product, "L1T") == 0)
        {
            if (strcmp (metadata->band[index].name, band_name) == 0)
            {
                input->thermal.filename =
                    strdup (metadata->band[index].file_name);

                input->thermal.toa_gain = metadata->band[index].toa_gain;
                input->thermal.toa_bias = metadata->band[index].toa_bias;

                input->thermal.size.s = metadata->band[index].nsamps;
                input->thermal.size.l = metadata->band[index].nlines;
                input->thermal.pixel_size[0] =
                    metadata->band[index].pixel_size[0];
                input->thermal.pixel_size[1] =
                    metadata->band[index].pixel_size[1];

                input->thermal.band_index = index;

                /* Only the one so break out */
                break;
            }
        }
    }

    /* Get the map projection coordinates */
    input->meta.ul_map_corner.x = metadata->global.proj_info.ul_corner[0];
    input->meta.ul_map_corner.y = metadata->global.proj_info.ul_corner[1];
    input->meta.ul_map_corner.is_fill = true;
    input->meta.lr_map_corner.x = metadata->global.proj_info.lr_corner[0];
    input->meta.lr_map_corner.y = metadata->global.proj_info.lr_corner[1];
    input->meta.lr_map_corner.is_fill = true;

    /* Get the geo graphic coordinates */
    input->meta.ul_geo_corner.lat = metadata->global.ul_corner[0];
    input->meta.ul_geo_corner.lon = metadata->global.ul_corner[1];
    input->meta.ul_geo_corner.is_fill = true;
    input->meta.lr_geo_corner.lat = metadata->global.lr_corner[0];
    input->meta.lr_geo_corner.lon = metadata->global.lr_corner[1];
    input->meta.lr_geo_corner.is_fill = true;

    /* Convert the acquisition date/time values */
    snprintf (acq_date, sizeof (acq_date), global->acquisition_date);
    snprintf (acq_time, sizeof (acq_time), global->scene_center_time);

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
