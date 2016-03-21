
#ifndef GET_ARGS_H
#define GET_ARGS_H


#include <stdbool.h>


int get_args
(
    int argc,           /* I: number of cmd-line args */
    char *argv[],       /* I: string of cmd-line args */
    char *xml_filename, /* I: address of input XML metadata filename  */
    bool *tape_6,       /* O: use the tape6 output */
    bool *verbose,      /* O: verbose flag */
    bool *debug         /* O: debug flag */
);


#endif /* GET_ARGS_H */
