
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <getopt.h>


#include "const.h"
#include "utilities.h"
#include "input.h"


/******************************************************************************
MODULE:  get_args

PURPOSE:  Gets the command-line arguments and validates that the required
arguments were specified.

RETURN VALUE:
Type = int
Value           Description
-----           -----------
FAILURE         Error getting the command-line arguments or a command-line
                argument and associated value were not specified
SUCCESS         No errors encountered

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/2/2013    Gail Schmidt     Original Development
3/15/2013   Song Guo         Changed to support Fmask
9/13/2013   Song Guo         Changed to use RETURN_ERROR
2/19/2014   Gail Schmidt     Modified to utilize the ESPA internal raw binary
                             file format

NOTES:
  1. Memory is allocated for the input and output files.  All of these should
     be character pointers set to NULL on input.  The caller is responsible
     for freeing the allocated memory upon successful return.
******************************************************************************/
int get_args
(
    int argc,                 /* I: number of cmd-line args */
    char *argv[],             /* I: string of cmd-line args */
    char **xml_infile,        /* I: address of input XML metadata filename  */
    char **dem_infile,        /* I: address of input DEM filename */
    char **emissivity_infile, /* I: address of input emissivity filename */
    bool *verbose             /* O: verbose flag */
)
{
    int c;                         /* current argument index */
    int option_index;              /* index for the command-line option */
    static int verbose_flag = 0;   /* verbose flag */
    char errmsg[MAX_STR_LEN];      /* error message */
    char FUNC_NAME[] = "get_args"; /* function name */
    static struct option long_options[] = {
        {"verbose", no_argument, &verbose_flag, 1},
        {"xml", required_argument, 0, 'i'},
        {"dem", required_argument, 0, 'd'},
        {"emi", required_argument, 0, 'e'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    /* Loop through all the cmd-line options */
    opterr = 0; /* turn off getopt_long error msgs as we'll print our own */
    while (1)
    {
        /* optstring in call to getopt_long is empty since we will only
           support the long options */
        c = getopt_long (argc, argv, "", long_options, &option_index);
        if (c == -1)
        {                       /* Out of cmd-line options */
            break;
        }

        switch (c)
        {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;

            case 'h':              /* help */
                usage ();
                return FAILURE;
                break;

            case 'i':              /* xml infile */
                *xml_infile = strdup (optarg);
                break;

            case 'd':              /* dem infile */
                *dem_infile = strdup (optarg);
                break;

            case 'e':              /* xml infile */
                *emissivity_infile = strdup (optarg);
                break;

            case '?':
            default:
                sprintf (errmsg, "Unknown option %s", argv[optind - 1]);
                usage ();
                RETURN_ERROR (errmsg, FUNC_NAME, FAILURE);
                break;
        }
    }

    /* Make sure the infile was specified */
    if (*xml_infile == NULL)
    {
        usage ();
        RETURN_ERROR ("XML input file is a required argument", FUNC_NAME,
                      FAILURE);
    }
    if (*dem_infile == NULL)
    {
        usage ();
        RETURN_ERROR ("DEM input file is a required argument", FUNC_NAME,
                      FAILURE);
    }
    if (*emissivity_infile == NULL)
    {
        usage ();
        RETURN_ERROR ("Emissivity input file is a required argument",
                      FUNC_NAME, FAILURE);
    }


    /* Check the verbose flag */
    if (verbose_flag)
        *verbose = true;
    else
        *verbose = false;

    if (*verbose)
    {
        printf ("XML_input_file = %s\n", *xml_infile);
        printf ("DEM_input_file = %s\n", *dem_infile);
        printf ("Emissivity_input_file = %s\n", *emissivity_infile);
    }

    return SUCCESS;
}
