#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

/* TODO:
 * If no data is given, i.e., if argc==3,
 * read from stdin
 */

void usage(char ** argv)
{
  printf("Usage:");
  printf("%s outfile format pd1 pd2 ...\n", argv[0]);
  printf(" where format is a c-type\n");
  printf("Example:\n");
  printf("%s R.double double 0.1 0.2\n", argv[0]);
  printf("%s W.uint8 uint8_t 0 1 1 0\n", argv[0]);
  return;
}

int main(int argc, char ** argv)
{
  if(argc<3)
  {
    usage(argv);
    return 1;
  }


  if(strcmp(argv[2], "uint8_t") == 0)
  {
  FILE * fout = fopen(argv[1], "wb");
  for(int kk = 3; kk<argc; kk++)
  {
    uint8_t val = atoi(argv[kk]);
    fwrite(&val, 1, sizeof(uint8_t), fout);
  }
  fclose(fout);
  return 0;
  }

   if(strcmp(argv[2], "double") == 0)
  {
  FILE * fout = fopen(argv[1], "wb");
  for(int kk = 3; kk<argc; kk++)
  {
    double val = atof(argv[kk]);
    fwrite(&val, 1, sizeof(double), fout);
  }
  fclose(fout);
  return 0;
  }


 printf("Does not recognized format: %s\n", argv[2]);
 printf("Supported formats: 'double', 'uint8_t'\n");

  return 0;
}
