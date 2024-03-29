/**
 * @file cf_util.c
 * @author Erik Wernersson
 * @date 2020-2023
 */

#include "cf_util.h"

int limit_mem(size_t max_bytes)
{
#ifdef __linux__
    if(max_bytes==0)
    {
        struct sysinfo info;
        if(sysinfo(&info))
        {
            fprintf(stderr, "sysinfo returned errno: %d\n", errno);
            errno = 0;
        }
        max_bytes = info.mem_unit*info.freeram;
    }
    struct rlimit limit;

    /* The limits are not really limits... */
    getrlimit(RLIMIT_DATA, &limit);
    limit.rlim_cur = max_bytes;
    setrlimit(RLIMIT_DATA, &limit);
    return EXIT_SUCCESS;
#else
    return EXIT_FAILURE;
#endif
}


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

double clockdiff(struct timespec* start,
                 struct timespec * finish)
{
    double elapsed = (finish->tv_sec - start->tv_sec);
    elapsed += (finish->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}

int64_t cf_file_size(const char * filename)
{
    FILE * fid = fopen(filename, "r");
    if(fid == NULL)
    {
        return -1;
    }
    fseek(fid, -0L, SEEK_END);
    size_t size = ftell(fid);
    fclose(fid);
    return (int64_t) size;
}
