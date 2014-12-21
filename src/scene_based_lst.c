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
    char *xml_name = NULL;              /* input XML filename */
    char *dem_name = NULL;              /* input DEM filename */
    char *emissivity_name = NULL;       /* input Emissivity filename */
    char directory[MAX_STR_LEN];        /* input/output data directory */
    char extension[MAX_STR_LEN];        /* input TOA file extension */
    Input_t *input = NULL;              /* input data and meta data */
    char scene_name[MAX_STR_LEN];       /* input data scene name */
    char command[MAX_STR_LEN];
    int status;                    /* return value from function call */
    //    Output_t *output = NULL; /* output structure and metadata */
    bool verbose;            /* verbose flag for printing messages */
    Espa_internal_meta_t xml_metadata;  /* XML metadata structure */
    float alb = 0.1;
    int i, k;
    int num_points;
    char **case_list = NULL;
    char **command_list = NULL;
    float **results = NULL;
    FILE *fd;
    int num_cases;
  
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

    /* Validate the input metadata file */
    if (validate_xml_file (xml_name) != SUCCESS)
    {  /* Error messages already written */
        LST_ERROR (errstr, "main");
    }

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
        printf ("DEBUG: Instrument: %d\n", input->meta.inst);
        printf ("DEBUG: Satellite: %d\n", input->meta.sat);
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

        printf("DEBUG: SUN AZIMUTH: %f\n", input->meta.sun_az);
        printf("DEBUG: SUN ZENITH: %f\n", input->meta.sun_zen);
        printf("DEBUG: Year, Month, Day, Hour, Minute, Second: %d, "
               "%d, %d, %d, %d,%f\n", input->meta.acq_date.year,
               input->meta.acq_date.month, input->meta.acq_date.day,
               input->meta.acq_date.hour, input->meta.acq_date.minute,
               input->meta.acq_date.second);
        printf("DEBUG: UL_MAP_CORNER: %f, %f\n", input->meta.ul_map_corner.x,
               input->meta.ul_map_corner.y);
        printf("DEBUG: LR_MAP_CORNER: %f, %f\n", input->meta.lr_map_corner.x,
               input->meta.lr_map_corner.y);
        printf("DEBUG: UL_GEO_CORNER: %f, %f\n", input->meta.ul_geo_corner.lat,
               input->meta.ul_geo_corner.lon);
        printf("DEBUG: LR_GEO_CORNER: %f, %f\n", input->meta.lr_geo_corner.lat,
               input->meta.lr_geo_corner.lon);
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
#endif
    /* Write out the intermediate values */
    fd = fopen("datetime.txt", "w"); 
    if (fd == NULL)
    {
        sprintf (errstr, "Opening file: datetime.txt\n");
        RETURN_ERROR (errstr, "scene_based_lst", FAILURE);
    }

    /* Write out the intermediate file */
    fprintf(fd, "%d,%d,%d,%d,%d,%f\n", input->meta.acq_date.year, 
            input->meta.acq_date.month, input->meta.acq_date.day, 
            input->meta.acq_date.hour, input->meta.acq_date.minute, 
            input->meta.acq_date.second);

    /* Close the intermediate file */
    status = fclose(fd);
    if ( status )
    {
        sprintf (errstr, "Closing file: datetime.txt\n");
        RETURN_ERROR (errstr, "scene_based_lst", FAILURE);
    }

    status = system("cp $BIN/narr_retrieval.bash .");
    if (status != SUCCESS)
    {
        sprintf (errstr, "cp $BIN/narr_retrieval.bash\n");
        RETURN_ERROR (errstr, "scene_based_lst", FAILURE);
    }
    
    status = system("./narr_retrieval.bash");
    if (status != SUCCESS)
    {
        sprintf (errstr, "./narr_retrieval.bash\n");
        RETURN_ERROR (errstr, "scene_based_lst", FAILURE);
    }    
    
    status = system("rm datetime.txt");
    if (status != SUCCESS)
    {
        sprintf (errstr, "rm datetime.txt\n");
        RETURN_ERROR (errstr, "scene_based_lst", FAILURE);
    }

    status = system("rm script*");
    if (status != SUCCESS)
    {
        sprintf (errstr, "rm script*\n");
        RETURN_ERROR (errstr, "scene_based_lst", FAILURE);
    }

    status = system("rm narr_retrieval.bash");
    if (status != SUCCESS)
    {
        sprintf (errstr, "rm narr_retrieval.bash\n");
        RETURN_ERROR (errstr, "scene_based_lst", FAILURE);
    }

    status = system("rm *_grb2txt");
    if (status != SUCCESS)
    {
        sprintf (errstr, "rm *_grb2txt\n");
        RETURN_ERROR (errstr, "scene_based_lst", FAILURE);
    }

    status = system("rm wgrib");
    if (status != SUCCESS)
    {
        sprintf (errstr, "\n");
        RETURN_ERROR (errstr, "scene_based_lst", FAILURE);
    }

    status = system("rm get_*.pl");
    if (status != SUCCESS)
    {
        sprintf (errstr, "rm get_*.pl\n");
        RETURN_ERROR (errstr, "scene_based_lst", FAILURE);
    }

    /* call first_files to generate tape5 file and commandList */
    status = first_files(input, &num_points, verbose);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Calling first_files\n");
        LST_ERROR (errstr, "scene_based_lst");
    }

    num_cases = num_points * NUM_ELEVATIONS * 3;
    case_list = (char **)allocate_2d_array(num_cases, MAX_STR_LEN, sizeof(char));  
    if (case_list == NULL)
    {
        sprintf (errstr, "Allocating case_list memory");
        LST_ERROR (errstr, "first_files");
    }

    command_list = (char **)allocate_2d_array(num_cases, MAX_STR_LEN, sizeof(char));  
    if (command_list == NULL)
    {
        sprintf (errstr, "Allocating command_list memory");
        LST_ERROR (errstr, "first_files");
    }

    /* Read case_list from caseList file */ 
    fd = fopen("caseList", "r"); 
    if (fd == NULL)
    {
        sprintf (errstr, "Opening file: caseList\n");
        LST_ERROR (errstr, "first_files");
    }

    /* Read in the caseList file */
    for (k = 0; k < num_cases; k++)
    {
        fscanf(fd, "%s", case_list[k]);
    }

    /* Close the caseList file */
    status = fclose(fd);
    if ( status )
    {
        sprintf (errstr, "Closing file: caseList\n");
        LST_ERROR (errstr, "first_files");
    }

    /* Read command_list from commandList file */
    fd = fopen("commandList", "r"); 
    if (fd == NULL)
    {
        sprintf (errstr, "Opening file: commandList\n");
        LST_ERROR (errstr, "first_files");
    }

    /* Read in the commandList file */
    for (k = 0; k < num_cases; k++)
    {
        fgets(command_list[k], MAX_STR_LEN, fd);
    }

    /* Close the commandList file */
    status = fclose(fd);
    if ( status )
    {
        sprintf (errstr, "Closing file: commandList\n");
        LST_ERROR (errstr, "first_files");
    }

    /* perform modtran runs by calling command_list */
    for (i = 0; i < num_cases; i++)
    {
        status = system(command_list[i]);
        if (status != SUCCESS)
        {
            sprintf (errstr, "executing command_list[i]");
            RETURN_ERROR (errstr, "scene_based_lst", FAILURE);
        }
    } 

    /* PARSING TAPE6 FILES: for each case in caseList (for each modtran run), 
       copy program to delete headers and parse wavelength and total radiance 
       from tape6 file*/
    status = system("cp $BIN/tape6parser.bash .");
    if (status != SUCCESS)
    {
        sprintf (errstr, "cp $BIN/tape6parser.bash\n");
        RETURN_ERROR (errstr, "scene_based_lst", FAILURE);
    }
    

    for (i = 0; i < num_cases; i++)
    {
        /* Just use $LST_DATA/elim2.sed directly instead of linking it */
        sprintf(command, "./tape6parser.bash %s", case_list[i]);
        status = system(command);
        if (status != SUCCESS)
        {
            sprintf (errstr, "./tape6parser.bash\n");
            RETURN_ERROR (errstr, "scene_based_lst", FAILURE);
        }
    }

    status = system("rm tape6parser.bash");
    if (status != SUCCESS)
    {
        sprintf (errstr, "rm tape6parser.bash\n");
        RETURN_ERROR (errstr, "scene_based_lst", FAILURE);
    }

    /* Free memory allocation */
    status = free_2d_array((void **)command_list);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: command_list\n");
        RETURN_ERROR (errstr, "scene_based_list", FAILURE);              
    }

    /* Allocate memory for results */
    results = (float **)allocate_2d_array(num_points * NUM_ELEVATIONS, 6,  
              sizeof(float)); 
    if (results == NULL)
    {
        sprintf (errstr, "Allocating results memory");
        LST_ERROR (errstr, "second_narr");
    }

    /* call second_narr to generate parameters for each height and NARR point */
    status = second_narr(input, num_points, alb, case_list, results, verbose);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Calling scene_based_list\n");
        LST_ERROR (errstr, "scene_based_lst");
    }

    /* Free memory allocation */
    status = free_2d_array((void **)case_list);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: current_case\n");
        RETURN_ERROR (errstr, "scene_based_list", FAILURE);              
    }

    /* call third_pixels_post to generate parameters for each Landsat pixel */
    status = third_pixels_post(input, num_points, dem_name, emissivity_name, results, 
             verbose);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Calling scene_based_list\n");
        LST_ERROR (errstr, "scene_based_lst");
    }

#if 0
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

    /* Close the input file and free the structure */
    CloseInput (input);
    FreeInput (input);

    free(xml_name);
    printf ("Processing complete.\n");
#endif

    /* Free memory allocation */
    status = free_2d_array((void **)results);
    if (status != SUCCESS)
    {
        sprintf (errstr, "Freeing memory: results\n");
        RETURN_ERROR (errstr, "scene_based_list", FAILURE);              
    }
#if 0
    /* Delete temporary file */
    status = system("rm newHead*");
    if (status != SUCCESS)
    {
        sprintf (errstr, "Deleting newHead* files\n");
        RETURN_ERROR (errstr, "scene_based_lst", FAILURE);
    }

    status = system("rm newTail*");
    if (status != SUCCESS)
    {
        sprintf (errstr, "Deleting newTail* files\n");
        RETURN_ERROR (errstr, "scene_based_lst", FAILURE);
    }

    status = system("rm tempLayers.txt");
    if (status != SUCCESS)
    {
        sprintf (errstr, "Deleting tempLayers file\n");
        RETURN_ERROR (errstr, "scene_based_lst", FAILURE);
    }

    /* Delete temporary directories */
    status = system("\rm -r HGT*");
    if (status != SUCCESS)
    {
        sprintf (errstr, "Deleting HGT* directories\n");
        RETURN_ERROR (errstr, "scene_based_lst", FAILURE);
    }

    status = system("\rm -r SHUM*");
    if (status != SUCCESS)
    {
        sprintf (errstr, "Deleting SHUM* directories\n");
        RETURN_ERROR (errstr, "scene_based_lst", FAILURE);
    }

    status = system("\rm -r TMP*");
    if (status != SUCCESS)
    {
        sprintf (errstr, "Deleting TMP* directories\n");
        RETURN_ERROR (errstr, "scene_based_lst", FAILURE);
    }

    status = system("\rm -r 4?.*_*");
    if (status != SUCCESS)
    {
        sprintf (errstr, "Deleting temporary directories\n");
        RETURN_ERROR (errstr, "scene_based_lst", FAILURE);
    }
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
    printf ("scene_base_lst \\n\n");
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
    printf ("\nExample: ./scene_based_lst --xml=LE70390032010263EDC00.xml "
            "--dem=17_30_DEM.tif "
            "--emi=AG100B.v003.-20.122.0001.bin --verbose\n");
    printf("Note: the scene_based_lst must run from the directory "
           "where the input data are located\n\n");
}
