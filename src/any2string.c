/**
 * @file any2string.c
 * @author Erik Wernersson
 * @date 2020-2023
 */

#include "any2string.h"

/* Would be nice:
 * -- Specify number of columns
 * -- Show summary
 * -- standard argument parsing
 */

static void usage(char ** argv)
{
    printf("Usage:\n");
    printf("%s format file\n", argv[0]);
    printf(" where format is a c-type\n");
    printf("Example:\n");
    printf("%s double R.double\n", argv[0]);
    printf("%s uint8_t L.uint8\n", argv[0]);
    return;
}

static FILE * openfile(const char * filename)
{
    FILE * fid = fopen(filename, "r");
    if(fid == NULL)
    {
        fprintf(stderr, "Unable to open %s for reading\n", filename);
        exit(EXIT_FAILURE);
    }
    return fid;
}
int any2string(int argc, char ** argv)
{
    if(argc<3)
    {
        usage(argv);
        return 1;
    }

    if(strcmp(argv[1], "double") == 0)
    {
        FILE * fin = openfile(argv[2]);
        double number = 0;
        while(fread(&number, sizeof(double), 1, fin) != 0)
        {
            printf("%f\n", number);
        }
        fclose(fin);
        return EXIT_SUCCESS;
    }

    if( (strcmp(argv[1], "uint8_t") == 0) | (strcmp(argv[1], "u8") == 0) )
    {
        FILE * fin = openfile(argv[2]);
        uint8_t byte = 0;
        while(fread(&byte, sizeof(uint8_t), 1, fin) != 0)
        {
            printf("%u\n", byte);
        }
        fclose(fin);
        return EXIT_SUCCESS;
    }

    if( (strcmp(argv[1], "uint16_t") == 0 ) | (strcmp(argv[1], "u16") == 0) )
    {
        FILE * fin = openfile(argv[2]);
        uint16_t token = 0;
        while(fread(&token, sizeof(uint16_t), 1, fin) != 0)
        {
            printf("%u\n", token);
        }
        fclose(fin);
        return EXIT_SUCCESS;
    }

    if( (strcmp(argv[1], "uint32_t") == 0) | (strcmp(argv[1], "u32") == 0) )
    {
        FILE * fin = openfile(argv[2]);
        uint32_t token = 0;
        while(fread(&token, sizeof(uint32_t), 1, fin) != 0)
        {
            printf("%u\n", token);
        }
        fclose(fin);
        return EXIT_SUCCESS;
    }

    printf("Does not recognized format: %s\n", argv[1]);
    printf("Supported formats: uint8_t, uint32_t, double\n");

    return EXIT_SUCCESS;
}
