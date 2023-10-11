#include "cf_util.h"

char * cf_timestr()
{
    size_t len = 128;
    char * str = calloc(len, 1);
    assert(str != NULL);

    time_t t = time(NULL);
    struct tm * tmp = localtime(&t);

    if (tmp == NULL) {
        perror("localtime");
        exit(EXIT_FAILURE);
    }

    if (strftime(str, len, "%F %H:%M:%S", tmp) == 0) {
        fprintf(stderr, "strftime returned 0");
        exit(EXIT_FAILURE);
    }
    return str;
}
