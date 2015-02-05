
#include <stdio.h>
#include <limits.h>


#include "const.h"
#include "utilities.h"


/******************************************************************************
MODULE:  read_narr_coordinates

PURPOSE: Creates directories and writes tape5 file, caseList, and commandList

RETURN: SUCCESS
        FAILURE
******************************************************************************/
int read_narr_coordinates
(
    char *lst_data_dir,
    int **grid_i,
    int **grid_j,
    float **lat,
    float **lon
)
{
    char FUNC_NAME[] = "read_narr_coordinates";
    int i, j;
    int count;
    char coord_file[PATH_MAX];
    FILE *fd = NULL;

    /* Setup the string to be used to open the coordinates file */
    count = snprintf (coord_file, sizeof (coord_file),
                      "%s/%s", lst_data_dir, "narr_coordinates.txt");
    if (count < 0 || count >= sizeof (coord_file))
    {
        RETURN_ERROR ("Failed initializing coord_file variable for"
                      " narr_coordinates.txt", FUNC_NAME, FAILURE);
    }

    /* Open the NARR coordinates file */
    fd = fopen (coord_file, "r");
    if (fd == NULL)
    {
        RETURN_ERROR ("Can't open narr_coordinates.txt file", FUNC_NAME,
                      FAILURE);
    }

    /* Read each value from the file */
    for (j = 0; j < NARR_COL; j++)
    {
        for (i = 0; i < NARR_ROW; i++)
        {
            /* File Format:
               Grid_I Grid_J Latitude Longitude
             */
            if (fscanf (fd, "%d %d %f %f",
                        &grid_i[i][j], &grid_j[i][j],
                        &lat[i][j], &lon[i][j]) == EOF)
            {
                RETURN_ERROR ("End of file (EOF) is met before"
                              " NARR_ROW * NARR_COL lines",
                              FUNC_NAME, FAILURE);
            }

            /* TODO - Should think about fixing the input file, so that this
                      confusing conversion is not needed.
                      When you read  the file data, it is as if you are
                      reading the values from the lower left to the upper
                      right as applied to the earth.  And the values being
                      read in start with an origin somewhere around the lower
                      left, hence the need for the following conversion.
               NOTE - If this is changed here, the else-where in the code will
                      break. */
            if ((lon[i][j] - 180.0) > MINSIGMA)
                lon[i][j] = 360.0 - lon[i][j];
            else
                lon[i][j] = -lon[i][j];
        }
    }

    if (fclose (fd) != SUCCESS)
    {
        RETURN_ERROR ("Closing file: narr_coordinates.txt\n",
                      FUNC_NAME, FAILURE);
    }

    return SUCCESS;
}
