#ifndef INTEGRATE_H
#define INTEGRATE_H

typedef struct
{
    int size;         /* size of arrays */
    double *array[2]; /* array pointers */
    int zsize;        /* size of z array */
    double *z;        /* interpolated values for integration */
    int isize;        /* size of integration point array */
    int *i;           /* integration points */
} integration_workspace;

int int_tabulated
(
    double *x,         /*I: Tabulated X-value data */
    double *f,         /*I: Tabulated F-value data */
    int nums,          /*I: Number of points */
    double *result_out /*O: Integrated result */
);

void free_integration_workspace();

#endif
