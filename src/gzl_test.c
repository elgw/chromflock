#include <stdio.h>
#include <stdlib.h>

#include "gzl.h"

/*
 * $ ./gzl_test R10F_contacts.pairs.gz  | md5sum
 * 59fbb98ca7f5326dc2e2f1b2f1106780  -
 * $ zcat R10F_contacts.pairs.gz | md5sum
 * 59fbb98ca7f5326dc2e2f1b2f1106780  -
 *
 */

int main(int argc, char ** argv)
{
    uint64_t buffer_size = 2048;
    if(argc < 2 || argc > 3)
    {
        printf("Usage: %s file.txt.gz [buffer_size]\n", argv[0]);
        return EXIT_FAILURE;
    }

    if(argc > 2)
    {
        buffer_size = atol(argv[2]);
    }

    gzl_state * gzl = gzl_open(argv[1], buffer_size);
    if(gzl == NULL)
    {
        fprintf(stderr, "Failed to open %s\n", argv[1]);
        return EXIT_FAILURE;
    }

    const char * str = NULL;
    int err = 0;
    while( (str = gzl_get_line(gzl, &err)) )
    {
        printf("%s\n", str);
    }
    int ret = EXIT_SUCCESS;
    if(err != 0)
    {
        fprintf(stderr, "gzl error code %d: %s\n",
                err, gzl_get_error_string(gzl));
        ret = EXIT_FAILURE;
    }
    gzl_destroy(gzl);
    return ret;
}
