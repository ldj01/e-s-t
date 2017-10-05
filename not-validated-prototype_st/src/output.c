
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <limits.h>
#include <time.h>


#include "espa_metadata.h"
#include "parse_metadata.h"
#include "write_metadata.h"
#include "envi_header.h"
#include "raw_binary_io.h"


#include "const.h"
#include "utilities.h"


#define MAX_DATE_LEN 28


/******************************************************************************
  NAME:  add_st_band_product

  PURPOSE:  Create a new envi output file including envi header and add the
            associated information to the XML metadata file.

  RETURN VALUE:  Type = int
      Value    Description
      -------  ---------------------------------------------------------------
      SUCCESS  No errors were encountered.
      ERROR    An error was encountered.
******************************************************************************/
int add_st_band_product
(
    char *xml_filename,
    char *reference_band_name,
    char *image_filename,
    char *product_name,
    char *band_name,
    char *short_name,
    char *long_name,
    char *data_units,
    int min_range,
    int max_range
)
{
    char FUNC_NAME[] = "add_st_band_product";

    int band_index = -1;
    int src_index = -1;
    char scene_name[PATH_MAX];
    char *tmp_char = NULL;
    Espa_internal_meta_t in_meta;
    Espa_internal_meta_t tmp_meta;
    Espa_band_meta_t *bmeta = NULL; /* pointer to the band metadata array
                                       within the output structure */
    time_t tp;                   /* time structure */
    struct tm *tm = NULL;        /* time structure for UTC time */
    char production_date[MAX_DATE_LEN+1]; /* current date/time for production */
    Envi_header_t envi_hdr;   /* output ENVI header information */
    char envi_file[PATH_MAX];

    /* Initialize the input metadata structure */
    init_metadata_struct (&in_meta);

    /* Parse the metadata file into our internal metadata structure; also
       allocates space as needed for various pointers in the global and band
       metadata */
    if (parse_metadata (xml_filename, &in_meta) != SUCCESS)
    {
        /* Error messages already written */
        return ERROR;
    }

    /* Find the representative band for metadata information */
    for (band_index = 0; band_index < in_meta.nbands; band_index++)
    {
        if (((strcmp (in_meta.band[band_index].product, "L1T") == 0)
             || (strcmp (in_meta.band[band_index].product, "L1G") == 0)
             || (strcmp (in_meta.band[band_index].product, "L1TP") == 0)
             || (strcmp (in_meta.band[band_index].product, "L1GT") == 0)
             || (strcmp (in_meta.band[band_index].product, "L1GS") == 0))
            && (strcmp (in_meta.band[band_index].name, reference_band_name)
                == 0))
        {
            /* this is the index we'll use for output band information */
            src_index = band_index;
            break;
        }
    }

    /* Figure out the scene name */
    strcpy (scene_name, in_meta.band[src_index].file_name);
    tmp_char = strchr (scene_name, '_');
    if (tmp_char != NULL)
        *tmp_char = '\0';

    /* Get the current date/time (UTC) for the production date of each band */
    if (time (&tp) == -1)
    {
        RETURN_ERROR ("unable to obtain current time", FUNC_NAME, ERROR);
    }

    tm = gmtime (&tp);
    if (tm == NULL)
    {
        RETURN_ERROR ("converting time to UTC", FUNC_NAME, ERROR);
    }

    if (strftime (production_date, MAX_DATE_LEN, "%Y-%m-%dT%H:%M:%SZ", tm) == 0)
    {
        RETURN_ERROR ("formatting the production date/time", FUNC_NAME, ERROR);
    }

    /* Gather all the band information from the representative band */

    /* Initialize the internal metadata for the output product. The global
       metadata won't be updated, however the band metadata will be updated
       and used later for appending to the original XML file. */
    init_metadata_struct (&tmp_meta);

    /* Allocate memory for the output band */
    if (allocate_band_metadata (&tmp_meta, 1) != SUCCESS)
        RETURN_ERROR("allocating band metadata", FUNC_NAME, ERROR);
    bmeta = tmp_meta.band;

    snprintf (bmeta[0].short_name, sizeof (bmeta[0].short_name),
              "%s", in_meta.global.product_id);
    bmeta[0].short_name[4] = '\0';
    strcat (bmeta[0].short_name, short_name);
    snprintf (bmeta[0].product, sizeof (bmeta[0].product), "%s",
              product_name);
    snprintf (bmeta[0].source, sizeof (bmeta[0].source), "level1");
    snprintf (bmeta[0].category, sizeof (bmeta[0].category), "image");
    bmeta[0].nlines = in_meta.band[src_index].nlines;
    bmeta[0].nsamps = in_meta.band[src_index].nsamps;
    bmeta[0].pixel_size[0] = in_meta.band[src_index].pixel_size[0];
    bmeta[0].pixel_size[1] = in_meta.band[src_index].pixel_size[1];
    snprintf (bmeta[0].pixel_units, sizeof (bmeta[0].pixel_units), "meters");
    snprintf (bmeta[0].app_version, sizeof (bmeta[0].app_version),
              "st_%s", ST_VERSION);
    snprintf (bmeta[0].production_date, sizeof (bmeta[0].production_date),
              "%s", production_date);
    bmeta[0].data_type = ESPA_FLOAT32;
    bmeta[0].fill_value = ST_NO_DATA_VALUE;
    bmeta[0].valid_range[0] = min_range;
    bmeta[0].valid_range[1] = max_range;
    snprintf (bmeta[0].name, sizeof (bmeta[0].name), "%s", band_name);
    snprintf (bmeta[0].long_name, sizeof (bmeta[0].long_name), "%s",
              long_name);
    snprintf (bmeta[0].data_units, sizeof (bmeta[0].data_units), "%s",
              data_units);
    snprintf (bmeta[0].file_name, sizeof (bmeta[0].file_name), "%s",
              image_filename);

    /* Create the ENVI header file for this band */
    if (create_envi_struct (&bmeta[0], &in_meta.global, &envi_hdr) != SUCCESS)
    {
        RETURN_ERROR ("Failed to create ENVI header structure.", FUNC_NAME,
                      ERROR);
    }

    /* Write the ENVI header */
    snprintf (envi_file, sizeof(envi_file), "%s", bmeta[0].file_name);
    tmp_char = strchr (envi_file, '.');
    if (tmp_char == NULL)
    {
        RETURN_ERROR ("Failed creating ENVI header filename", FUNC_NAME,
                      ERROR);
    }

    sprintf (tmp_char, ".hdr");
    if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
    {
        RETURN_ERROR ("Failed writing ENVI header file", FUNC_NAME, ERROR);
    }

    /* Append the ST band to the XML file */
    if (append_metadata (1, bmeta, xml_filename)
        != SUCCESS)
    {
        RETURN_ERROR ("Appending spectral index bands to XML file",
                       FUNC_NAME, ERROR);
    }

    free_metadata (&in_meta);
    free_metadata (&tmp_meta);

    return SUCCESS;
}
