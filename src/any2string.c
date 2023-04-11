#include "any2string.h"

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

int any2string(int argc, char ** argv)
{
    if(argc<3)
    {
        usage(argv);
        return 1;
    }


    if(strcmp(argv[1], "double") == 0)
    {
        FILE * fin = fopen(argv[2], "r");
        double number = 0;
        while(fread(&number, sizeof(double), 1, fin) != 0)
        {
            printf("%f\n", number);
        }
        fclose(fin);
        return 0;
    }

    if(strcmp(argv[1], "uint8_t") == 0)
    {
        FILE * fin = fopen(argv[2], "r");
        uint8_t byte = 0;
        while(fread(&byte, sizeof(uint8_t), 1, fin) != 0)
        {
            printf("%u\n", byte);
        }
        fclose(fin);
        return EXIT_SUCCESS;
    }

    printf("Does not recognized format: %s\n", argv[1]);
    printf("Supported formats: 'uint8_t'\n");

    return EXIT_SUCCESS;
}
