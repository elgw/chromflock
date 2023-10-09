#include "string2any.h"

/* TODO:
 * If no data is given, i.e., if argc==3,
 * read from stdin
 */


static void usage(char ** argv)
{
    printf("Usage:\n");
    printf("%s outfile format pd1 pd2 ...\n", argv[0]);
    printf(" where format is a c-type\n");
    printf("Example:\n");
    printf("%s R.double double 0.1 0.2\n", argv[0]);
    printf("%s W.uint8 uint8_t 0 1 1 0\n", argv[0]);
    return;
}

int string2any(int argc, char ** argv)
{
    if(argc<3)
    {
        usage(argv);
        return EXIT_FAILURE;
    }

    if(strcmp(argv[2], "uint8_t") == 0)
    {
        FILE * fout = fopen(argv[1], "wb");
        if(fout == NULL)
        {
            fprintf(stderr, "Unable to open %s for writing\n", argv[1]);
            exit(EXIT_FAILURE);
        }
        for(int kk = 3; kk<argc; kk++)
        {
            uint8_t val = atoi(argv[kk]);
            fwrite(&val, 1, sizeof(uint8_t), fout);
        }
        fclose(fout);
        return EXIT_SUCCESS;
    }

    if(strcmp(argv[2], "uint32_t") == 0)
    {
        FILE * fout = fopen(argv[1], "wb");
        if(fout == NULL)
        {
            fprintf(stderr, "Unable to open %s for writing\n", argv[1]);
            exit(EXIT_FAILURE);
        }
        for(int kk = 3; kk<argc; kk++)
        {
            uint32_t val = atoi(argv[kk]);
            fwrite(&val, 1, sizeof(uint32_t), fout);
        }
        fclose(fout);
        return EXIT_SUCCESS;
    }

    if(strcmp(argv[2], "double") == 0)
    {
        FILE * fout = fopen(argv[1], "wb");
        if(fout == NULL)
        {
            fprintf(stderr, "Unable to open %s for writing\n", argv[1]);
            exit(EXIT_FAILURE);
        }
        for(int kk = 3; kk<argc; kk++)
        {
            double val = atof(argv[kk]);
            fwrite(&val, 1, sizeof(double), fout);
        }
        fclose(fout);
        return EXIT_SUCCESS;
    }


    fprintf(stderr, "Does not recognized format: %s\n", argv[2]);
    fprintf(stderr, "Supported formats: 'double', 'uint8_t' and 'uint32_t'\n");
    fprintf(stderr, "\n");
    usage(argv);

    return EXIT_FAILURE;
}
