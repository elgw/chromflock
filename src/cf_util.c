/**
 * @file cf_util.c
 * @author Erik Wernersson
 * @date 2020-2023
 */


#include <ctype.h>

#include "cf_util.h"
#include "npio.h"

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

const char* cf_YES_NO(int v)
{
    if(v == 1)
    { return "YES"; }
    return "NO";
}

void bpos_print(FILE * fid, bpos * P)
{
    fprintf(fid, "#=%u (x=%f, y=%f, z=%f)\n", P->bead_id,
            P->x, P->y, P->z);
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
                               float * X32, double * X64,
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
    if(X32)
    {
        for( ; ll < nbead; ll++)
        {
            int read = getline(&line, &line_len, fid);
            if(read == -1){ goto parsing_error; }
            char *ptr = strtok(line, delim);
            if(ptr == NULL){ goto parsing_error; }
            X32[3*ll] = atof(ptr);
            ptr = strtok(NULL, delim);
            if(ptr == NULL){ goto parsing_error; }
            X32[3*ll+1] = atof(ptr);
            ptr = strtok(NULL, delim);
            if(ptr == NULL){ goto parsing_error; }
            X32[3*ll+2] = atof(ptr);
        }
    }
    if(X64)
    {
        for( ; ll < nbead; ll++)
        {
            int read = getline(&line, &line_len, fid);
            if(read == -1){ goto parsing_error; }
            char *ptr = strtok(line, delim);
            if(ptr == NULL){ goto parsing_error; }
            X64[3*ll] = atof(ptr);
            ptr = strtok(NULL, delim);
            if(ptr == NULL){ goto parsing_error; }
            X64[3*ll+1] = atof(ptr);
            ptr = strtok(NULL, delim);
            if(ptr == NULL){ goto parsing_error; }
            X64[3*ll+2] = atof(ptr);
        }
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

int
load_bead_coordinates_from_npy(const char * fname,
                               float * X32, double * X64,
                               const i64 nbead)
{
    npio_t * npy = npio_load(fname);
    if(npy == NULL)
    {
        return -1;
    }
    if(npy->ndim != 2)
    {
        goto fail_data;
    }

    if(npy->shape[1] != 4)
    {
        goto fail_data;
    }

    if(npy->dtype != NPIO_F32)
    {
        goto fail_data;
    }
    if(X64 != NULL)
    {
        const float * C = (const float *) npy->data;
        for(i64 kk = 0; kk < nbead; kk++)
        {
            for(i64 ll =0; ll < 3; ll++)
            {
                X64[3*kk+ll] = C[4*kk+ll];
            }
        }
    }
    if(X32 != NULL)
    {
        const float * C = (const float *) npy->data;
        for(i64 kk = 0; kk < nbead; kk++)
        {
            for(i64 ll =0; ll < 3; ll++)
            {
                X32[3*kk+ll] = C[4*kk+ll];
            }
        }
    }
    npio_free(npy);
    return 0;

 fail_data:
    fprintf(stderr, "%s contains invalid data\n", fname);
    npio_print(stderr, npy);
    npio_free(npy);
    return -1;
}

bpos *
load_bead_apos_from_npy(const char * fname,
                        int * nconstraint)
{
    npio_t * npy = npio_load(fname);
    if(npy == NULL)
    {
        return NULL;
    }
    if(npy->ndim != 2)
    {
        goto fail_data;
    }

    int nrow = (int) npy->shape[0];
    int ncol = (int) npy->shape[1];
    if(ncol != 4)
    {
        goto fail_data;
    }

    if(npy->dtype == NPIO_F32)
    {

        *nconstraint = nrow;
        bpos * apos = calloc(nrow, sizeof(bpos));
        if(apos == NULL)
        {
            printf("Memory allocation failure while reading %s\n", fname);
            exit(EXIT_FAILURE);
        }

        const float * C = (const float *) npy->data;
        for(i64 kk = 0; kk < nrow; kk++)
        {
            apos[kk].bead_id = (u32) C[4*kk];
            apos[kk].x = C[4*kk + 1];
            apos[kk].y = C[4*kk + 2];
            apos[kk].z = C[4*kk + 3];
            float r = pow(apos[kk].x, 2) + pow(apos[kk].y, 2) + pow(apos[kk].z, 2);
            if(r > 1)
            {
                bpos_print(stdout, apos+kk);
            }
        }

        npio_free(npy);
        return apos;
    }
    fprintf(stderr, "Invalid data format, expected F32\n");

 fail_data:
    fprintf(stderr, "%s contains invalid data\n", fname);
    npio_print(stderr, npy);
    npio_free(npy);
    return NULL;
}

u32 *
load_bead_contacts_from_npy(const char * fname, int * _ncont)
{
    npio_t * npy = npio_load(fname);
    if(npy == NULL)
    {
        printf("Could not open %s as an npy file\n", fname);
        return NULL;
    }
    if(npy->nel == 0)
    {
        npio_free(npy);
        *_ncont = 0;
        return malloc(0); // ok on linux 6.7
    }
    if(npy->ndim != 2)
    {
        goto fail_data;
    }

    if(npy->shape[1] != 2)
    {
        goto fail_data;
    }
    int ncont = npy->shape[0];
    *_ncont = ncont;

    if(npy->dtype == NPIO_U32)
    {
        u32 * C = calloc(2*ncont, sizeof(u32));
        assert(C != NULL);
        memcpy(C, npy->data, 2*ncont*sizeof(u32));
        npio_free(npy);
        return C;
    }
    if(npy ->dtype == NPIO_I64)
    {
        u32 * C = calloc(2*ncont, sizeof(u32));
        assert(C != NULL);
        i64 * _C = (i64*) npy->data;
        for(i64 kk = 0; kk < 2*ncont; kk++)
        {
            C[kk] = _C[kk];
        }
        npio_free(npy);
        return C;
    }

 fail_data:
    fprintf(stderr, "%s contains invalid data\n", fname);
    npio_print(stderr, npy);
    npio_free(npy);
    return NULL;
}

u8 *
load_bead_labels_from_npy(const char * fname, int * _nbead)
{
    npio_t * npy = npio_load(fname);
    if(npy == NULL)
    {
        printf("Could not open %s as an npy file\n", fname);
        return NULL;
    }

    int nbead = npy->nel;
    if(nbead == 0)
    {
        fprintf(stderr, "Empty array");
        goto fail_data;
    }
    *_nbead = nbead;

    if(npy->dtype == NPIO_U8)
    {
        u8 * L = calloc(nbead, 1);
        if(L == NULL)
        {
            exit(EXIT_FAILURE);
        }
        memcpy(L, npy->data, nbead);
        npio_free(npy);
        return L;
    }

    if(npy->dtype == NPIO_I64)
    {
        u8 * L = calloc(nbead, 1);
        if(L == NULL)
        {
            exit(EXIT_FAILURE);
        }
        i64 * in_L = (i64*) npy->data;
        for(i64 kk = 0; kk < nbead; kk++)
        {
            L[kk] = (u8) in_L[kk];
        }
        npio_free(npy);
        return L;
    }

 fail_data:
    fprintf(stderr, "%s contains invalid data\n", fname);
    npio_print(stderr, npy);
    npio_free(npy);
    return NULL;
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
                              const mflock_geometry_type geometry,
                              const elli * ell)
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
        if(ell == NULL)
        {
            radius = norm3d(X+3*kk);
        } else {
            radius = elli_getScale(ell, X+3*kk);
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


int write_bead_coordinates_to_npy(const char * fname,
                                  const double * X,
                                  const i64 nbead,
                                  const mflock_geometry_type geometry,
                                  const elli * ell)
{
    float * C = calloc(4*nbead, sizeof(float));
    if(C == NULL)
    {
        fprintf(stderr,
                "write_bead_coordinates_to_npy: "
                "memory allocation failure\n");
        return -1;
    }

    i64 n_outside = 0;

    for(i64 kk = 0; kk<nbead; kk++)
    {
        for(i64 ll = 0; ll < 3; ll++)
        {
            C[4*kk + ll] = X[3*kk + ll];
        }
        double radius;
        if( (ell != NULL) & geometry == MFLOCK_ELLIPSOID) {
            radius = elli_getScale(ell, X+3*kk);
        } else {
            radius = norm3d(X+3*kk);
        }

        C[4*kk + 3] = radius;

        switch(geometry)
        {
        case MFLOCK_SPHERE:
        case MFLOCK_ELLIPSOID:
            if(radius > 1.0)
            {
                n_outside++;
            }
            break;
        case MFLOCK_BOX:
            {
                int inside = 1;
                for(int bb = 0; bb < 3; bb++)
                {
                    if(C[bb] > 1.0 | C[bb] < -1.0) {
                        inside = 0;
                    }
                }

                if(inside == 0)
                {
                    n_outside++;
                }
            }
        }
    }

    if(n_outside > 0)
    {
        printf("Warning %ld beads are outside of the domain\n", n_outside);
    }

    int shape[2] = {nbead, 4};

    i64 status = npio_write(fname, 2,
                            shape,
                            (const void *) C,
                            NPIO_F32, NPIO_F32);
    free(C);
    if(status == -1)
    {
        fprintf(stderr,
                "write_bead_coordinates_to_npy: "
                "npio_write failed when writing to %s\n", fname);
        return -1;
    }

    return 0;
}


int npy_extension(const char * name)
{
    assert(name != NULL);
    if(name == NULL)
    {
        return 0;
    }
    size_t len = strlen(name);
    if(len < 5)
    {
        return 0;
    }
    if(toupper(name[len-1]) != 'Y')
    {
        return 0;
    }
    if(toupper(name[len-2]) != 'P')
    {
        return 0;
    }
    if(toupper(name[len-3]) != 'N')
    {
        return 0;
    }
    if(name[len-4] != '.')
    {
        return 0;
    }
    return 1;
}

int
write_bead_labels(const char * fname, const u8 * labels, i64 nbin)
{
    if(npy_extension(fname))
    {
        return write_bead_labels_to_npy(fname, labels, nbin);
    } else {
        return write_bead_labels_to_u8(fname, labels, nbin);
    }
}

int
write_bead_labels_to_u8(const char * fname,
                        const u8 * labels,
                        i64 nbin)
{

    FILE * fid = fopen(fname, "wb");
    if(fid == NULL)
    {
        fprintf(stderr, "Error opening %s\n", fname);
        return -1;
    }
    size_t nwritten = fwrite(labels, sizeof(u8), nbin, fid);
    fclose(fid);
    if((i64) nwritten != nbin)
    {
        fprintf(stderr, "Error writing to %s\n", fname);
        return -1;
    }
    return 0;
}


int
write_bead_labels_to_npy(const char * fname,
                         const u8 * labels,
                         i64 nbin)
{
    assert(nbin >= 0);
    assert(fname != NULL);
    assert(labels != NULL);
    int shape[] = {nbin};
    i64 nwritten =  npio_write(fname,
                               1, shape,
                               (const void *) labels,
                               NPIO_U8, NPIO_U8);
    if(nwritten > 0)
    {
        return 0;
    }
    return -1;
}

int
write_u32(const char * fname, u32 * data,
          i64 M, i64 N)
{
    if(npy_extension(fname))
    {
        return write_u32_to_npy(fname, data, M, N);
    } else {
        return write_u32_to_raw(fname, data, M, N);
    }
}

int write_u32_to_raw(const char * fname,
                     u32 * data,
                     i64 M, i64 N)
{
    FILE * fid = fopen(fname, "wb");
    if(fid == NULL)
    {
        return -1;
    }
    size_t nwritten = fwrite(data, sizeof(u32), M*N, fid);
    if(nwritten != (size_t) M*N)
    {
        fclose(fid);
        return -1;
    }
    fclose(fid);
    return 0;
}

int write_u32_to_npy(const char * fname,
                     u32 * data,
                     i64 M, i64 N)
{
    assert(fname != NULL);
    assert(data != NULL);
    assert(M>0);
    assert(N>0);

    int shape[] = {N, M}; // flip dimensions
    i64 nwritten = npio_write(fname,
                              2, shape,
                              (const void *) data,
                              NPIO_U32, NPIO_U32);
    if(nwritten > 0)
    {
        return 0;
    } else {
        return -1;
    }
}

u8 * load_cmap(const char * fname)
{
    npio_t * npy = npio_load(fname);
    if(npy == NULL)
    {
        return NULL;
    }
    if(npy->ndim != 2)
    {
        printf("The color map is not a 2D array\n");
        npio_free(npy);
        return NULL;
    }

    if(npy->shape[1] != 3)
    {
        printf("The color map has the wrong shape, expected 256x3, got %dx%d\n",
               npy->shape[1], npy->shape[0]);
        npio_free(npy);
        return NULL;
    }

    if(npy->shape[0] != 256)
    {
        printf("The color map has the wrong shape, expected 256x3, got %dx%d\n",
               npy->shape[1], npy->shape[0]);
        npio_free(npy);
        return NULL;
    }

    if(npy->dtype != NPIO_U8)
    {
        printf("The color map has the wrong data type, should be uint8\n");
        npio_free(npy);
        return NULL;
    }

    uint8_t * cmap = npy->data;
    npy->data = NULL;
    npio_free(npy);
    return cmap;
}
