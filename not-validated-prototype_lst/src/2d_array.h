
#ifndef MISC_2D_ARRAY_H
#define MISC_2D_ARRAY_H

/*****************************************************************************
  DESCRIPTION:  Allocates and provides access to 2 dimensional arrays.

  HISTORY:  This has derived and modified from the original implementation
            found in Landsat 8 IAS source.
*****************************************************************************/

#include <stdio.h>


void **allocate_2d_array
(
    int rows,          /* I: Number of rows for the 2D array */
    int columns,       /* I: Number of columns for the 2D array */
    size_t member_size /* I: Size of the 2D array element */
);


int free_2d_array
(
    void **array_ptr   /* I: Pointer returned by the alloc routine */
);


#endif
