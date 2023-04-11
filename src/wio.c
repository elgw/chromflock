#include "wio.h"

/* read compressed */
static void*
wio_read_z(char * fileName, size_t * nel);

/* read uncompressed */
static void*
wio_read_u(char * fileName, size_t * nel);

/* write compressed */
static int
wio_write_z(char * fileName, size_t nel, void * data);

/* write uncompressed */
static int
wio_write_u(char * fileName, size_t nel, void * data);

void * wio_read(char * fileName, size_t * nel)
{
    int fnl = strlen(fileName);

    // If very short file name, it does not end with .gz
    if(fnl < 3)
    {
        return wio_read_u(fileName, nel);
    }

    // If file name ends with '.gz'
//  printf("%d %c%c%c\n", fnl, fileName[fnl-3], fileName[fnl-2], fileName[fnl-1]);
    if(fileName[fnl-3] == '.')
    {
        if(fileName[fnl-2] == 'g')
        {
            if(fileName[fnl-1] == 'z')
            {
                return wio_read_z(fileName, nel);
            }
        }
    }

    return wio_read_u(fileName, nel);
}

int wio_write(char * fileName, size_t nel, void * data)
{

    int fnl = strlen(fileName);

    // If very short file name, it does not end with .gz
    if(fnl < 3)
    {
        return wio_write_u(fileName, nel, data);
    }

    // If file name ends with '.gz'
//  printf("%d %c%c%c\n", fnl, fileName[fnl-3], fileName[fnl-2], fileName[fnl-1]);
    if(fileName[fnl-3] == '.')
    {
        if(fileName[fnl-2] == 'g')
        {
            if(fileName[fnl-1] == 'z')
            {
                return wio_write_z(fileName, nel, data);
            }
        }
    }

    return wio_write_u(fileName, nel, data);
}

int wio_write_z(char * fileName, size_t nel, void * data)
{

    if(nel > INT_MAX)
    {
        printf("gzwrite can't write a file of size %zu, max size is %u\n", nel, INT_MAX);
        exit(1);
    }

    printf("wio_write_z(%s, %zu, %ptr)\n", fileName, nel, data);
    gzFile zf = gzopen(fileName, "wb1");
    if(zf == Z_NULL)
    {
        printf("gzopen failed for file %s\n", fileName);
        exit(1);
    }
    int status = gzwrite(zf, data, nel);

    if(status <= 0) // 0 bytes written
    {
        int z_errnum = 0;
        printf("gzwrite error %s\n", gzerror(zf, &z_errnum));
        exit(1);
    }
    gzclose(zf);

    return 0;
}

int wio_write_u(char * fileName, size_t nel, void * restrict data)
{
    FILE * f = fopen(fileName, "wb1");
    fwrite((void *) data, sizeof(uint8_t), nel, f);
    fclose(f);
    return 0;
}


void * wio_read_z(char * fileName, size_t * nel)
{
    /* TODO: Can we use inflate() and deflate() from zlib to
     * read/write files larger than 4 GB ?. Even better is of course
     * to store the W matrices as sparse data */

    gzFile zf = gzopen(fileName, "rb");
    // ftell...


    int fd = open(fileName, O_RDONLY);
    if (fd == -1) {
        printf("wio_read_z: Error opening %s\n", fileName);
        return NULL;
    }

    // Last four bytes gives the file size mod 2**32
    // so this isn't perfectly correct
    lseek(fd, -4L, SEEK_END);
    uint32_t size = 0;
    size_t read_byte = read(fd, &size, 4); // Read 20 bytes
    if(read_byte != 4)
    {
        printf("Couldn't read the last 4 bytes\n");
    }
    //printf("File size mod 2*32 = %u\n", size);
    close(fd);

    nel[0] = size;

    uint8_t * w = malloc(nel[0]*sizeof(uint8_t));
    if (w == NULL) {
        fprintf(stderr, "Fatal: failed to allocate %zu bytes.\n", nel[0]);
        abort();
    }
    gzread(zf, w, nel[0]*sizeof(uint8_t));
    gzclose(zf);
    return w;
}

void * wio_read_u(char * fileName, size_t * nel)
{
    FILE * fd = fopen(fileName, "r");
    if (fd == NULL) {
        printf("wio_read_u: Error opening file\n");
        return NULL;
    }

    // Last four bytes gives the file size mod 2**32
    // so this isn't perfectly correct
    fseek(fd, -0L, SEEK_END);
    size_t size = ftell(fd);
    nel[0] = size;
    //printf("File size = %zu\n", size);
    rewind(fd);

    uint8_t * w = malloc(size*sizeof(uint8_t));

    if (w == NULL) {
        fprintf(stderr, "Fatal: failed to allocate %zu bytes.\n", size);
        abort();
    }
    size_t read_byte = fread(w, 1, size, fd);
    if(read_byte != size)
    {
        printf("%s/%s/%d Failed reading %s\n",
               __FILE__, __FUNCTION__, __LINE__, fileName);
        free(w);
        return NULL;
    }

    fclose(fd);
    return w;
}



int wio_ut(int argc, char ** argv)
{

    /* Create a uint8_t matrix
     * write it both compressed and non-compressed
     * read back both copies
     * and compare with the original
     */

    if(argc != 1)
    {
        printf("%s should not be called with any arguments\n", argv[0]);
    }

    char wuFile[] = "w.test.uint8";
    char wzFile[] = "w.test.uint8.gz";

    size_t N = 3000;
    uint8_t * w = malloc(N*N*sizeof(uint8_t));
    size_t nset = 0;
    for(size_t kk  = 0; kk<N*N; kk++)
    {
        if( rand() > .9*RAND_MAX)
        {
            w[kk] = 1;
            nset++;
        } else
        {
            w[kk] = 0;
        }
    }
    printf("%zu / %zu set to 1\n", nset, N*N);

    wio_write(wuFile, N*N, w); // write uncompressed
    wio_write(wzFile, N*N, w); // write compressed

    size_t nelz = 0;
    uint8_t * wz = wio_read(wzFile, &nelz);

    size_t nelu = 0;
    uint8_t * wu = wio_read(wuFile, &nelu);

    printf("nelz: %zu nelu: %zu\n", nelz, nelu);
    size_t diffu = 0;
    size_t diffz = 0;
    for(size_t kk = 0; kk<nelz; kk++)
    {
        if(w[kk] != wz[kk])
            diffz++;
        if(w[kk] != wu[kk])
            diffu++;
    }
    printf("diffu: %zu diffz: %zu\n", diffu, diffz);

    assert(diffu == 0);
    assert(diffz == 0);

    free(wz);
    free(wu);
    free(w);

    return EXIT_SUCCESS;
}
