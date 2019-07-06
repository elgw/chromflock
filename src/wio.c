#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdint.h>
#include <assert.h>

// read compressed
void * wio_read_z(char * fileName, size_t * nel);
// read uncompressed
void * wio_read_u(char * fileName, size_t * nel);

// write compressed
int wio_write_z(char * fileName, size_t nel, void * data);
// write uncompressed
int wio_write_u(char * fileName, size_t nel, void * data);

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

  gzFile zf = gzopen(fileName, "wb");
  gzwrite(zf, data, nel);
  gzclose(zf);

  return 0;
}

int wio_write_u(char * fileName, size_t nel, void * restrict data)
{
  FILE * f = fopen(fileName, "wb");
  fwrite((void *) data, sizeof(uint8_t), nel, f);
  fclose(f);
  return 0;
}


void * wio_read_z(char * fileName, size_t * nel)
{
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
  gzread(zf, w, nel[0]*sizeof(uint8_t));
  gzclose(zf);
  return w;
}

void * wio_read_u(char * fileName, size_t * nel)
{
  int fd = open(fileName, O_RDONLY);
  if (fd == -1) {
    printf("Error opening file\n");
    return NULL;
  }

  // Last four bytes gives the file size mod 2**32
  // so this isn't perfectly correct
  size_t size = lseek(fd, -0L, SEEK_END);
  nel[0] = size;
 // printf("File size = %zu\n", size);
  lseek(fd, 0, SEEK_SET); // rewind

  uint8_t * w = malloc(size*sizeof(uint8_t));
  size_t read_byte = read(fd, w, size);
  if(read_byte != size)
  {
    printf("Failed reading the file\n");
  }

  close(fd);
  return w;

}


#ifdef WIO_UT
int main(int argc, char ** argv)
{

  /* Create a uint8_t matrix
   * write it both compressed and non-compressed
   * read back both copies
   * and compare with the original
   */

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

  return 0;
}
#endif
