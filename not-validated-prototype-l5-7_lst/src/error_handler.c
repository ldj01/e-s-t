#include "swe.h"

/******************************************************************************
MODULE:  error_handler

PURPOSE:  Prints the error/warning message.

RETURN VALUE:
Type = None

PROJECT:  Land Satellites Data System Science Research and Development (LSRD)
at the USGS EROS

HISTORY:
Date        Programmer       Reason
--------    ---------------  -------------------------------------
1/2/2012    Gail Schmidt     Original Development
4/26/2013   Song Guo         For surface water extent use

NOTES:
******************************************************************************/
void error_handler
(
    bool error_flag,  /* I: true for errors, false for warnings */
    char *module,     /* I: calling module name */
    char *errmsg      /* I: error message to be printed, without ending EOL */
)
{
    if (error_flag)
        printf ("Error: %s : %s\n\n", module, errmsg);
    else
        printf ("Warning: %s : %s\n", module, errmsg);
}
