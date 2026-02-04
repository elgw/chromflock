/**
 * @file cf_util.c
 * @author Erik Wernersson
 * @date 2020-2023
 */

#include "cf_util.h"

typedef int64_t i64;

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

int
load_bead_coordinates_from_csv(const char * fname,
                               double * X,
                               const i64 nbead)
{

    FILE * fid = fopen(fname, "r");
    if(fid == NULL)
    {
        fprintf(stderr, "Can't open %s for reading\n", fname);
        return -1;
    }

    size_t line_len = 1024;
    char * line = calloc(line_len, sizeof(char));
    if(line == NULL)
    {
        fprintf(stderr, "Memory allocation error\n");
        return -1;
    }

    char delim[] = ",";
    i64 ll = 0;
    for( ; ll < nbead; ll++)
    {
        int read = getline(&line, &line_len, fid);
        if(read == -1){ goto parsing_error; }
        char *ptr = strtok(line, delim);
        if(ptr == NULL){ goto parsing_error; }
        X[3*ll] = atof(ptr);
        ptr = strtok(NULL, delim);
        if(ptr == NULL){ goto parsing_error; }
        X[3*ll+1] = atof(ptr);
        ptr = strtok(NULL, delim);
        if(ptr == NULL){ goto parsing_error; }
        X[3*ll+2] = atof(ptr);
    }

    fclose(fid);
    free(line);
    return 0;

 parsing_error:
    fclose(fid);
    fprintf(stderr, "Failed to read line %zu from %s\n", ll+1, fname);
    free(line);
    return -1;

}

static double norm3d(const double * restrict X)
{
    double n = 0;
    for(size_t kk = 0; kk<3; kk++)
        n+=pow(X[kk], 2);
    return sqrt(n);
}


int
write_bead_coordinates_to_csv(const char * fname,
                              const double * X,
                              const i64 nbead,
                              const elli * geometry)
{
    FILE * fid = fopen(fname, "w");
    if(fid == NULL)
    {
        fprintf(stderr, "Unable to open %s for writing\n", fname);
        return -1;
    }
    for(i64 kk = 0; kk<nbead; kk++)
    {
        double radius;
        if(geometry == NULL)
        {
            radius = norm3d(X+3*kk);
        } else {
            radius = elli_getScale(geometry, X+3*kk);
        }

        int nwritten = fprintf(fid, "%f, %f, %f, %f\n",
                               X[3*kk], X[3*kk+1], X[3*kk+2],
                               radius);
        // TODO: Unless we check the length of the string to
        // write, we don't know if all bytes were written
        if(nwritten <= 0)
        {
            goto fail_write;
        }

    }
    fclose(fid);
    return 0;

 fail_write:
    fprintf(stderr, "An error occurred while writing to %s\n", fname);
    fclose(fid);
    return -1;
}
