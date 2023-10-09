#include "contact_pairs_io.h"

/* Check if the file ends with .gz */
static int isgzfile(const char * file)
{
    int fnl = strlen(file);
    if(fnl > 2)
    {
        if(file[fnl-3] == '.')
        {
            if(file[fnl-2] == 'g')
            {
                if(file[fnl-1] == 'z')
                {
                    return 1;
                }
            }
        }
    }
    return 0;
}

static uint32_t * contact_pairs_read_raw(const char * filename,
                                         uint64_t * nPairs)
{
    FILE * fid = fopen(filename, "r");
    if(fid == NULL)
    {
        return NULL;
    }
    fseek(fid, -0L, SEEK_END);
    size_t size = ftell(fid);
    rewind(fid);

    if( size % 2*sizeof(uint32_t) != 0)
    {
        fprintf(stderr, "%s has an unexpected size and might be corrupt\n",
                filename);
        fclose(fid);
        return NULL;
    }

    uint32_t * P = malloc(size);
    if(P == NULL)
    {
        fclose(fid);
        return NULL;
    }

    size_t nread = fread(P, sizeof(uint32_t),
                         size/sizeof(uint32_t), fid);
    if(nread != size/sizeof(uint32_t) )
    {
        fprintf(stderr, "Unable to read %zu numbers from %s, only got %zu\n",
                size/sizeof(uint32_t), filename, nread);
        free(P);
        fclose(fid);
        return NULL;
    }

    fclose(fid);
    *nPairs = nread/2;
    return P;
}

static size_t get_gz_size(const char * file)
{
    int fd = open(file, O_RDONLY);
    if (fd == -1) {
        fprintf(stderr, "%s %d: Error opening %s\n",
                __FILE__, __LINE__, file);
        return 0;
    }

    // Last four bytes gives the file size mod 2**32
    // so this isn't perfectly correct
    lseek(fd, -4L, SEEK_END);
    uint32_t size = 0;
    size_t read_byte = read(fd, &size, 4); // Read 20 bytes
    if(read_byte != 4)
    {
        fprintf(stderr, "Couldn't read the last 4 bytes from %s\n",
                file);
        close(fd);
        return 0;
    }
    //printf("File size mod 2*32 = %u\n", size);
    close(fd);

    return size;
}

static uint32_t * contact_pairs_read_gz(const char * file,
                                     uint64_t * nPairs)
{
    size_t uncompressed_size = get_gz_size(file);

    if(uncompressed_size == 0)
    {
        fprintf(stderr, "Unable to determine the uncompressed size of %s\n",
                file);
        return NULL;
    }

    uint32_t * CP = calloc(uncompressed_size/sizeof(uint32_t), sizeof(uint32_t));
    if(CP == NULL)
    {
        fprintf(stderr, "%s %d -- memory allocation failed\n",
                __FILE__, __LINE__);
        return NULL;
    }

    gzFile zf = gzopen(file, "rb");
    if(zf == Z_NULL)
    {
        fprintf(stderr, "gzopen failed for file %s\n", file);;
        return NULL;
    }

    int status = gzread(zf, CP, uncompressed_size);
    if(status < 1)
    {
        fprintf(stderr, "%s %d\n", __FILE__, __LINE__);
        fprintf(stderr, "gzerror %d\n", status);
        gzclose(zf);
        free(CP);
        return NULL;
    }
    gzclose(zf);
    *nPairs = uncompressed_size / sizeof(uint32_t) / 2;
    return CP;
}

uint32_t * contact_pairs_read(const char * file,
                                     uint64_t * nPairs)
{
    uint32_t * CP = NULL;
    if(isgzfile(file))
    {
        CP = contact_pairs_read_gz(file, nPairs);
    } else {
        CP = contact_pairs_read_raw(file, nPairs);
    }

    #if 0
    if(CP != NULL)
    {
        uint32_t max = CP[0];
        for(size_t kk = 0; kk< *nPairs*2; kk++)
        {
            CP[kk] > max ? max = CP[kk] : 0;
        }
        printf("%lu Contact pairs indicate at least %u beads\n", *nPairs, max);
    }
#endif

    return CP;
}

static int contact_pairs_write_raw(const char * file,
                                   const uint32_t * P,
                                   uint64_t nPairs)
{
    FILE * fid = fopen(file, "w");
    if(fid == NULL)
    {
        fprintf(stderr, "Unable to open %s for writing\n", file);
        return EXIT_FAILURE;
    }

    size_t nwritten = fwrite(P, sizeof(uint32_t), nPairs*2, fid);
    if(nwritten != 2*nPairs)
    {
        fprintf(stderr, "%s %d\n", __FILE__, __LINE__);
        fprintf(stderr, "Unable to write to %s\n", file);
        fprintf(stderr, "nwritten = %zu\n", nwritten);
        fclose(fid);
        return EXIT_FAILURE;
    }

    fclose(fid);
    return EXIT_SUCCESS;
}

static int contact_pairs_write_gz(const char * file,
                                  const uint32_t * P,
                                  uint64_t nPairs)
{
    size_t nbytes = nPairs*2*sizeof(uint32_t);
    if(nbytes > INT_MAX)
    {
        printf("gzwrite can't write a file of size %zu, max size is %u\n",
               nbytes, INT_MAX);
        return EXIT_FAILURE;
    }
    gzFile zf = gzopen(file, "wb1");
    if(zf == NULL)
    {
        fprintf(stderr, "Unable to open %s for writing\n", file);
        return EXIT_FAILURE;
    }
    int status = gzwrite(zf, P, nbytes);
    if(status < 1)
    {
        fprintf(stderr, "gzwrite failed for %s\n", file);
        gzclose(zf);
        return EXIT_FAILURE;
    }
    gzclose(zf);
    return EXIT_SUCCESS;
}

int contact_pairs_write(const char * file,
                        const uint32_t * P,
                        uint64_t nPairs)
{
    if(isgzfile(file))
    {
        return contact_pairs_write_gz(file, P, nPairs);
    } else
    {
        return contact_pairs_write_raw(file, P, nPairs);
    }
}

int contact_pairs_io_ut(int argc, char ** argv)
{
    if(argc != 1)
    {
        printf("%s should not be called with any arguments\n", argv[0]);
    }

    /* Simple test. Allocate, write, read, compare to allocated */

    size_t nPairs = 123;
    uint32_t * CP = calloc(nPairs*2, sizeof(uint32_t));
    for(size_t kk = 0; kk<nPairs*2; kk++)
    {
        CP[kk] = rand() % UINT32_MAX;
    }

    printf(" -> Writing and reading uncompressed data\n");
    char tempfile[] = "contact_pairs_io_ut_XXXXXX.u32";
    int file = mkstemps(tempfile, 4);
    if(file == -1)
    {
        fprintf(stderr, "%s %d\n", __FILE__, __LINE__);
        fprintf(stderr, "Unable to create temp file %s\n", tempfile);
        exit(EXIT_FAILURE);
    }
    close(file);
    if(contact_pairs_write_raw(tempfile, CP, nPairs))
    {
        fprintf(stderr, "contact_pairs_write_raw filed\n");
        exit(EXIT_FAILURE);
    }

    uint64_t nPairs2;
    uint32_t * CP_raw = contact_pairs_read(tempfile, &nPairs2);
    if(nPairs != nPairs2)
    {
        fprintf(stderr, "Wrote %lu pairs but read %lu pairs\n", nPairs, nPairs2);
        exit(EXIT_FAILURE);
    }
    size_t nDiff = 0;
    for(size_t kk = 0; kk<nPairs*2; kk++)
    {
        if(CP[kk] != CP_raw[kk])
        {
            nDiff++;
        }
    }
    if(nDiff > 0)
    {
        fprintf(stderr, "contact_pairs_write / contact_pairs_read failed\n");
        exit(EXIT_FAILURE);
    }

    printf(" -> Writing and reading compressed data\n");
    char tempfilegz[] = "contact_pairs_io_ut_XXXXXX.u32.gz";
    file = mkstemps(tempfilegz, 7);
    if(file == -1)
    {
        fprintf(stderr, "%s %d\n", __FILE__, __LINE__);
        fprintf(stderr, "Unable to create temp file %s\n", tempfilegz);
        exit(EXIT_FAILURE);
    }
    contact_pairs_write_gz(tempfilegz, CP, nPairs);
    uint64_t nPairs3;
    uint32_t * CP_gz = contact_pairs_read(tempfilegz, &nPairs3);
    if(nPairs != nPairs3)
    {
        fprintf(stderr, "Wrote %lu pairs but read %lu pairs\n", nPairs, nPairs3);
        exit(EXIT_FAILURE);
    }
    nDiff = 0;
    for(size_t kk = 0; kk<nPairs*2; kk++)
    {
        if(CP[kk] != CP_gz[kk])
        {
            nDiff++;
        }
    }
    if(nDiff > 0)
    {
        fprintf(stderr, "contact_pairs_write_gz / contact_pairs_read_gz failed\n");
        exit(EXIT_FAILURE);
    }

    free(CP);
    free(CP_raw);
    free(CP_gz);
    if(remove(tempfile))
    {
        fprintf(stderr, "Failed to remove %s\n", tempfile);
    }
    if(remove(tempfilegz))
    {
        fprintf(stderr, "Failed to remove %s\n", tempfilegz);
    }
    printf("All tests passed and temporary files cleaned up\n");


    return EXIT_SUCCESS;
}
