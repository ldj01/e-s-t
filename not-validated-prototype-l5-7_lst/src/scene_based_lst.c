#include <string.h>
#include <stdarg.h>
#include <time.h>
#include "input.h"
#include "output.h"
#include "scene_based_lst.h"
#include "2d_array.h"
#include "const.h"

/******************************************************************************
METHOD:  scene_based_lst

PURPOSE:  the main routine for scene based LST (Land Surface Temperature) in C

RETURN VALUE:
Type = int
Value           Description
-----           -----------
ERROR           An error occurred during processing of the scene_based_lst
SUCCESS         Processing was successful

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
3/15/2013   Song Guo         Original Development

NOTES: type ./scene_based_lst --help for information to run the code
******************************************************************************/
int main (int argc, char *argv[])
{
    char errstr[MAX_STR_LEN];           /* error string */
    char *cptr = NULL;                  /* pointer to the file extension */
    char *xml_name = NULL;              /* input XML filename */
    char *dem_name = NULL;              /* input DEM filename */
    char *emissivity_name = NULL;       /* input Emissivity filename */
    char directory[MAX_STR_LEN];        /* input/output data directory */
    char extension[MAX_STR_LEN];        /* input TOA file extension */
    Input_t *input = NULL;              /* input data and meta data */
    char envi_file[MAX_STR_LEN];        /* output ENVI file name */
    char scene_name[MAX_STR_LEN];       /* input data scene name */
    unsigned char **pixel_mask;    /* pixel mask */
    int status;                    /* return value from function call */
    Output_t *output = NULL; /* output structure and metadata */
    bool verbose;            /* verbose flag for printing messages */
    float sun_azi_temp = 0.0;/* Keep the original sun azimuth angle */
    Espa_internal_meta_t xml_metadata;  /* XML metadata structure */
    Envi_header_t envi_hdr;   /* output ENVI header information */
  
    time_t now;
    time(&now);
    printf("scene_based_lst start_time=%s\n",ctime(&now));

    /* Read the command-line arguments, including the name of the input
       Landsat TOA reflectance product and the DEM */
    status = get_args (argc, argv, &xml_name, &dem_name, &emissivity_name, 
                       &verbose);
    if (status != SUCCESS)
    { 
        sprintf (errstr, "calling get_args");
        LST_ERROR (errstr, "main");
    }
#if 0
    /* Validate the input metadata file */
    if (validate_xml_file (xml_name) != SUCCESS)
    {  /* Error messages already written */
        LST_ERROR (errstr, "main");
    }
#endif
    /* Initialize the metadata structure */
    init_metadata_struct (&xml_metadata);

    /* Parse the metadata file into our internal metadata structure; also
       allocates space as needed for various pointers in the global and band
       metadata */
    if (parse_metadata (xml_name, &xml_metadata) != SUCCESS)
    {  /* Error messages already written */
        LST_ERROR (errstr, "main");
    }

    /* Split the filename to obtain the directory, scene name, and extension */
    split_filename(xml_name, directory, scene_name, extension);
    if (verbose)
        printf("directory, scene_name, extension=%s,%s,%s\n", 
            directory, scene_name, extension);

    /* Open input file, read metadata, and set up buffers */
    input = OpenInput(&xml_metadata);
    if (input == NULL)
    {
        sprintf (errstr, "opening the TOA and brightness temp files in: %s",
                 xml_name);
        LST_ERROR (errstr, "main");
    }

    if (verbose)
    {
        /* Print some info to show how the input metadata works */
        printf ("DEBUG: Number of input thermal bands: %d\n", input->nband_th);
        printf ("DEBUG: Number of input lines: %d\n", input->size_th.l);
        printf ("DEBUG: Number of input samples: %d\n", input->size_th.s);
        printf ("DEBUG: ACQUISITION_DATE.DOY is %d\n", input->meta.acq_date.doy);
        printf ("DEBUG: Fill value is %d\n", input->meta.fill);
        printf ("DEBUG: Thermal Band -->\n");
        printf ("DEBUG:   therm_satu_value_ref: %d\n", 
                   input->meta.therm_satu_value_ref);
        printf ("DEBUG:   therm_satu_value_max: %d\n", 
                   input->meta.therm_satu_value_max);
        printf ("DEBUG:   therm_gain: %f, therm_bias: %f\n", 
                input->meta.gain_th, input->meta.bias_th);

        printf("DEBUG: SUN AZIMUTH is %f\n", input->meta.sun_az);
        printf("DEBUG: SUN ZENITH is %f\n", input->meta.sun_zen);
    }

#if 0
    /* If the scene is an ascending polar scene (flipped upside down), then
       the solar azimuth needs to be adjusted by 180 degrees.  The scene in
       this case would be north down and the solar azimuth is based on north
       being up clock-wise direction. Flip the south to be up will not change 
       the actual sun location, with the below relations, the solar azimuth
       angle will need add in 180.0 for correct sun location */
    if (input->meta.ul_corner.is_fill &&
        input->meta.lr_corner.is_fill &&
        (input->meta.ul_corner.lat - input->meta.lr_corner.lat) < MINSIGMA)
    {
        /* Keep the original solar azimuth angle */
        sun_azi_temp = input->meta.sun_az;
        input->meta.sun_az += 180.0;
        if ((input->meta.sun_az - 360.0) > MINSIGMA)
            input->meta.sun_az -= 360.0;
        if (verbose)
            printf ("  Polar or ascending scene.  Readjusting solar azimuth by "
                "180 degrees.\n  New value: %f degrees\n", input->meta.sun_az);
    }

    /* Retrieve the NARR data */
    status = system("narr_retrieval.bash input->meta.acq_date.year input->meta.acq_date.month input->meta.acq_date.day input->meta.acq_date.hour input->meta.acq_date.minute input->meta.acq_date.second");
    if (status != SUCCESS)
    {
        RETURN_ERROR (errstr, "pcloud", FAILURE);
    }    


    /* Reassign solar azimuth angle for output purpose if south up north 
       down scene is involved */
    if (input->meta.ul_corner.is_fill &&
        input->meta.lr_corner.is_fill &&
        (input->meta.ul_corner.lat - input->meta.lr_corner.lat) < MINSIGMA)
        input->meta.sun_az = sun_azi_temp;
   
    /* Open the output file */
    output = OpenOutput (&xml_metadata, input);
    if (output == NULL)
    {   /* error message already printed */
        sprintf(errstr, "Opening output file");
        LST_ERROR (errstr, "main");
    }

    if (!PutOutput(output, pixel_mask))
    {
        sprintf (errstr, "Writing output LST in HDF files\n");
        LST_ERROR (errstr, "main");
    }

    /* Close the output file */
    if (!CloseOutput (output))
    {
        sprintf (errstr, "closing output file");
        LST_ERROR(errstr, "main");
    }

    /* Create the ENVI header file this band */
    if (create_envi_struct (&output->metadata.band[0], &xml_metadata.global,
        &envi_hdr) != SUCCESS)
    {
        sprintf (errstr, "Creating ENVI header structure.");
        LST_ERROR(errstr, "main");
    }

    /* Write the ENVI header */
    strcpy (envi_file, output->metadata.band[0].file_name);
    cptr = strchr (envi_file, '.');
    if (cptr == NULL)
    {
        sprintf (errstr, "error in ENVI header filename");
        LST_ERROR(errstr, "main");
    }

    strcpy (cptr, ".hdr");
    if (write_envi_hdr (envi_file, &envi_hdr) != SUCCESS)
    {
        sprintf (errstr, "Writing ENVI header file.");
        LST_ERROR(errstr, "main");
    }

    /* Append the LST band to the XML file */
    if (append_metadata (output->nband, output->metadata.band, xml_name)
        != SUCCESS)
    {
        sprintf (errstr, "Appending spectral index bands to XML file.");
        LST_ERROR(errstr, "main");
    }
  
    /* Free the structure */
    if (!FreeOutput (output))
    {
        sprintf (errstr, "freeing output file structure");
        LST_ERROR(errstr, "main");
    }

    /* Free the metadata structure */
    free_metadata (&xml_metadata);

    /* Free the pixel mask */
    status = free_2d_array((void **)pixel_mask);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing mask memory");
        LST_ERROR(errstr, "main");
    }

    /* Close the input file and free the structure */
    CloseInput (input);
    FreeInput (input);

    free(xml_name);
    printf ("Processing complete.\n");
#endif
    time(&now);
    printf("scene_based_lst end_time=%s\n",ctime(&now));
    return (SUCCESS);
}

/******************************************************************************
MODULE:  usage

PURPOSE:  Prints the usage information for this application.

RETURN VALUE:
Type = None

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
3/15/2013   Song Guo         Original Development
8/15/2013   Song Guo         Modified to use TOA reflectance file 
                             as input instead of metadata file
2/19/2014   Gail Schmidt     Modified to utilize the ESPA internal raw binary
                             file format

NOTES: 
******************************************************************************/
void usage ()
{
    printf ("LST identify the cloud, shadow, snow, water and clear pixels using "
            "the input Landsat scene (top of atmosphere (TOA) reflection and "
            "brightness temperature (BT) for band 6) output from LEDAPS\n\n");
    printf ("usage: ./scene_based_lst "
            "--xml=input_xml_filename "
            "--dem=input_dem_filename "
            "--emi=input_emissivity_filename "
            "[--verbose]\n");

    printf ("\nwhere the following parameters are required:\n");
    printf ("    -xml: name of the input XML file\n");
    printf ("\nwhere the following parameters are optional:\n");
    printf ("    -verbose: should intermediate messages be printed? (default "
            "is false)\n");
    printf ("\n./scene_based_lst --help will print the usage statement\n");
    printf ("\nExample: ./cfmask --xml=LE70390032010263EDC00.xml "
            "--dem 17_30_DEM.tif "
            "--emi AG100B.v003.-20.122.0001.bin --verbose\n");
}
