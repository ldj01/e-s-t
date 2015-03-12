
#include <string.h>
#include <stdarg.h>
#include <unistd.h>
#include <time.h>


#include "const.h"
#include "2d_array.h"
#include "utilities.h"
#include "get_args.h"
#include "input.h"
#include "output.h"
#include "build_points.h"
#include "build_modtran_input.h"
#include "calculate_point_atmospheric_parameters.h"
#include "calculate_pixel_atmospheric_parameters.h"


/* These are for compile time debugging logic.
   Set them to 0 to turn them off.
   They need to be set to 1 for production/standard processing. */
#define RUN_MODTRAN 1
#define EXTRACT_TAPE6_RESULTS 1


/******************************************************************************
METHOD:  lst

PURPOSE:  The main routine for scene based LST (Land Surface Temperature).

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           An error occurred during processing of the scene_based_lst
SUCCESS         Processing was successful

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
          at the USGS EROS
******************************************************************************/
int
main (int argc, char *argv[])
{
    char FUNC_NAME[] = "main";

    Espa_internal_meta_t xml_metadata;  /* XML metadata structure */

    char msg_str[MAX_STR_LEN];
    char xml_filename[PATH_MAX];        /* input XML filename */
    char dem_filename[PATH_MAX];        /* input DEM filename */
    char emissivity_filename[PATH_MAX]; /* input Emissivity filename */
    char command[PATH_MAX];

    Input_t *input = NULL;          /* input data and meta data */
    //    Output_t *output = NULL; /* output structure and metadata */

    bool use_tape6;             /* Use the tape6 output */
    bool verbose;               /* verbose flag for printing messages */
    bool debug;                 /* debug flag for debug output */

    int modtran_run;

    double alb = 0.1;
    double **modtran_results = NULL;

    char *tmp_env = NULL;

    REANALYSIS_POINTS points;

    time_t now;

    /* Display the starting time of the application */
    time (&now);
    snprintf (msg_str, sizeof(msg_str),
              "LST start_time [%s]", ctime (&now));
    LOG_MESSAGE (msg_str, FUNC_NAME);

    /* Read the command-line arguments, including the name of the input
       Landsat TOA reflectance product and the DEM */
    if (get_args (argc, argv, xml_filename, dem_filename, emissivity_filename,
                  &use_tape6, &verbose, &debug) != SUCCESS)
    {
        RETURN_ERROR ("calling get_args", FUNC_NAME, EXIT_FAILURE);
    }

    /* Verify the existence of required environment variables */
    /* Grab the environment path to the LST_DATA_DIR */
    tmp_env = getenv ("LST_DATA_DIR");
    if (tmp_env == NULL)
    {
        RETURN_ERROR ("LST_DATA_DIR environment variable is not set",
                      FUNC_NAME, EXIT_FAILURE);
    }

    /* Validate the input metadata file */
    if (validate_xml_file (xml_filename) != SUCCESS)
    {
        /* Error messages already written */
        return EXIT_FAILURE;
    }

    /* Initialize the metadata structure */
    init_metadata_struct (&xml_metadata);

    /* Parse the metadata file into our internal metadata structure; also
       allocates space as needed for various pointers in the global and band
       metadata */
    if (parse_metadata (xml_filename, &xml_metadata) != SUCCESS)
    {
        /* Error messages already written */
        return EXIT_FAILURE;
    }

    /* Open input file, read metadata, and set up buffers */
    input = OpenInput (&xml_metadata);
    if (input == NULL)
    {
        RETURN_ERROR ("opening input files", FUNC_NAME, EXIT_FAILURE);
    }

    if (verbose)
    {
        /* Print some info to show how the input metadata works */
        printf ("Satellite: %d\n", input->meta.satellite);
        printf ("Instrument: %d\n", input->meta.instrument);

        printf ("Number of input lines: %d\n", input->thermal.size.l);
        printf ("Number of input samples: %d\n", input->thermal.size.s);

        printf ("Fill value is %d\n", input->thermal.fill_value);

        printf ("Thermal Band -->\n");
        printf ("  therm_gain: %f\n  therm_bias: %f\n",
                input->thermal.toa_gain, input->thermal.toa_bias);

        printf ("Year, Month, Day, Hour, Minute, Second:"
                " %d, %d, %d, %d, %d, %f\n",
                input->meta.acq_date.year, input->meta.acq_date.month,
                input->meta.acq_date.day, input->meta.acq_date.hour,
                input->meta.acq_date.minute, input->meta.acq_date.second);
        printf ("ACQUISITION_DATE.DOY is %d\n",
                input->meta.acq_date.doy);

        printf ("UL_MAP_CORNER: %f, %f\n", input->meta.ul_map_corner.x,
                input->meta.ul_map_corner.y);
        printf ("LR_MAP_CORNER: %f, %f\n", input->meta.lr_map_corner.x,
                input->meta.lr_map_corner.y);
        printf ("UL_GEO_CORNER: %f, %f\n",
                input->meta.ul_geo_corner.lat, input->meta.ul_geo_corner.lon);
        printf ("LR_GEO_CORNER: %f, %f\n",
                input->meta.lr_geo_corner.lat, input->meta.lr_geo_corner.lon);
    }

    /* Build the points that will be used */
    if (build_points (input, &points) != SUCCESS)
    {
        RETURN_ERROR ("Building POINTS input\n", FUNC_NAME, EXIT_FAILURE);
    }

    if (verbose)
    {
        printf ("Number of Points: %d\n", points.num_points);
    }

    /* Call build_modtran_input to generate the tape5 file input and
       the MODTRAN commands for each point and height */
    if (build_modtran_input (input, &points, verbose, debug)
        != SUCCESS)
    {
        RETURN_ERROR ("Building MODTRAN input\n", FUNC_NAME, EXIT_FAILURE);
    }

    /* Perform MODTRAN runs by calling each command */
    for (modtran_run = 0; modtran_run < points.num_modtran_runs; modtran_run++)
    {
        snprintf (msg_str, sizeof(msg_str),
                  "Executing MODTRAN [%s]",
                   points.modtran_runs[modtran_run].command);
        LOG_MESSAGE (msg_str, FUNC_NAME);

#if RUN_MODTRAN
        if (system (points.modtran_runs[modtran_run].command) != SUCCESS)
        {
            RETURN_ERROR ("Error executing MODTRAN", FUNC_NAME,
                          EXIT_FAILURE);
        }
#endif
    }

    /* PARSING MODTRAN RESULTS:
       for each case in caseList (for each modtran run),
       parse wavelength and total radiance from tape6 file into parsed */
    for (modtran_run = 0; modtran_run < points.num_modtran_runs; modtran_run++)
    {
        if (use_tape6)
        {
            /* Use modtran generated tape6 output */
            snprintf (command, sizeof (command),
                      "lst_extract_modtran_results.py"
                      " --tape6"
                      " --input-path %s"
                      " --output-path %s",
                      points.modtran_runs[modtran_run].path,
                      points.modtran_runs[modtran_run].path);
        }
        else
        {
            /* Use modtran generated pltout.asc output */
            snprintf (command, sizeof (command),
                      "lst_extract_modtran_results.py"
                      " --pltout"
                      " --input-path %s"
                      " --output-path %s",
                      points.modtran_runs[modtran_run].path,
                      points.modtran_runs[modtran_run].path);
        }

        snprintf (msg_str, sizeof(msg_str), "Executing [%s]", command);
        LOG_MESSAGE (msg_str, FUNC_NAME);

#if EXTRACT_TAPE6_RESULTS
        if (system (command) != SUCCESS)
        {
            RETURN_ERROR ("Failed executing lst_extract_tape6_results.py",
                          FUNC_NAME, EXIT_FAILURE);
        }
#endif
    }

    /* Allocate memory for MODTRAN results */
    modtran_results =
        (double **) allocate_2d_array (points.num_points * NUM_ELEVATIONS,
                                       MGPE_NUM_ELEMENTS, sizeof (double));
    if (modtran_results == NULL)
    {
        RETURN_ERROR ("Allocating MODTRAN results memory", FUNC_NAME,
                      EXIT_FAILURE);
    }

    /* Generate parameters for each height and NARR point */
    if (calculate_point_atmospheric_parameters (input, &points, alb,
                                                modtran_results, verbose)
        != SUCCESS)
    {
        RETURN_ERROR ("Calculating point atmospheric parameters\n",
                      FUNC_NAME, EXIT_FAILURE);
    }

    /* Generate parameters for each Landsat pixel */
    if (calculate_pixel_atmospheric_parameters (input, &points,
                                                xml_filename,
                                                dem_filename,
                                                emissivity_filename,
                                                modtran_results, verbose)
        != SUCCESS)
    {
        RETURN_ERROR ("Calculating per/pixel atmospheric parameters\n",
                      FUNC_NAME, EXIT_FAILURE);
    }

    /* Free memory allocation */
    free_points_memory (&points);

#if NOT_TESTED
    /* Open the output file */
    output = OpenOutput (&xml_metadata, input);
    if (output == NULL)
    {                           /* error message already printed */
        RETURN_ERROR ("Opening output file", FUNC_NAME, EXIT_FAILURE);
    }

    if (!PutOutput (output, pixel_mask))
    {
        RETURN_ERROR ("Writing output LST in HDF files\n", FUNC_NAME,
                      EXIT_FAILURE);
    }

    /* Close the output file */
    if (!CloseOutput (output))
    {
        RETURN_ERROR ("closing output file", FUNC_NAME, EXIT_FAILURE);
    }

    /* Create the ENVI header file this band */
    if (create_envi_struct (&output->metadata.band[0], &xml_metadata.global,
                            &envi_hdr) != SUCCESS)
    {
        RETURN_ERROR ("Creating ENVI header structure.", FUNC_NAME,
                      EXIT_FAILURE);
    }

    /* Write the ENVI header */
    strcpy (envi_file, output->metadata.band[0].file_name);
    cptr = strchr (envi_file, '.');
    if (cptr == NULL)
    {
        RETURN_ERROR ("error in ENVI header filename", FUNC_NAME,
                      EXIT_FAILURE);
    }

    strcpy (cptr, ".hdr");
    if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
    {
        RETURN_ERROR ("Writing ENVI header file.", FUNC_NAME, EXIT_FAILURE);
    }

    /* Append the LST band to the XML file */
    if (append_metadata (output->nband, output->metadata.band, xml_filename)
        != SUCCESS)
    {
        RETURN_ERROR ("Appending spectral index bands to XML file.",
                      FUNC_NAME, EXIT_FAILURE);
    }

    /* Free the structure */
    if (!FreeOutput (output))
    {
        RETURN_ERROR ("freeing output file structure", FUNC_NAME,
                      EXIT_FAILURE);
    }
#endif

    /* Free the metadata structure */
    free_metadata (&xml_metadata);

    /* Close the input file and free the structure */
    CloseInput (input);
    FreeInput (input);

    /* Free memory allocations */
    if (free_2d_array ((void **) modtran_results) != SUCCESS)
    {
        RETURN_ERROR ("Freeing memory: MODTRAN results\n", FUNC_NAME,
                      EXIT_FAILURE);
    }

    if (!debug)
    {
        /* Delete temporary file */
        if (unlink ("atmospheric_parameters.txt") != SUCCESS)
        {
            RETURN_ERROR ("Deleting atmospheric_parameters.txt files\n",
                          FUNC_NAME, EXIT_FAILURE);
        }

        if (unlink ("base_head.txt") != SUCCESS)
        {
            RETURN_ERROR ("Deleting baseHead.txt files\n", FUNC_NAME,
                          EXIT_FAILURE);
        }

        if (unlink ("new_tail.txt") != SUCCESS)
        {
            RETURN_ERROR ("Deleting newTail.txt files\n", FUNC_NAME,
                          EXIT_FAILURE);
        }

        if (unlink ("temp_layers.txt") != SUCCESS)
        {
            RETURN_ERROR ("Deleting tempLayers.txt file\n", FUNC_NAME,
                          EXIT_FAILURE);
        }

        if (unlink ("used_points.txt") != SUCCESS)
        {
            RETURN_ERROR ("Deleting used_points.txt file\n", FUNC_NAME,
                          EXIT_FAILURE);
        }
    }

    time (&now);
    snprintf (msg_str, sizeof(msg_str),
              "scene_based_lst end_time=%s\n", ctime (&now));
    LOG_MESSAGE (msg_str, FUNC_NAME);

    return EXIT_SUCCESS;
}
