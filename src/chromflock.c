#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>

#include "cf_version.h"
#include "chromflock_init.h"
#include "cc2cpm.h"
#include "string2any.h"
#include "any2string.h"
#include "sprite2cmap.h"

static int show_version(void)
{
    printf("chromflock v.%s\n", cf_version);
    return EXIT_SUCCESS;
}

static int usage(void)
{
    printf("usage:\n\t"
           "chromflock <command>\n");
    printf("\n");
    printf("Available commands:\n");
    printf("help\n\t"
           "Show this help message\n");
    printf("init\n\t"
           "Initialize this directory for chromflock structures\n");
    printf("hic2cpm\n\t"
           "Generate contact a contact probability matrix from Hi-C data\n");
    printf("sprite2cpm\n\t"
           "convert sprite date (.cluster files) to contact probability maps\n");
    printf("any2string\n\t"
           "convert raw uint8 and double to text\n");
    printf("string2any\n\t"
           "create small test data sets from strings to raw\n");
    printf("version\n\t"
           "show version information\n");
    printf("Each command has a separate help section\n");
    printf("\n");
    printf("Web page: https://www.github.com/elgw/chromflock\n");
    return EXIT_SUCCESS;
}


int main(int argc, char ** argv)
{

    if(argc == 1)
    {
        usage();
        return EXIT_SUCCESS;
    }

    const char * command = argv[1];

    if(!strcmp(command, "help")
       || !strcmp(command, "--help")
       || !strcmp(command, "-h"))
    {
        return usage();
    }

    if(!strcmp(command, "init"))
    {
        return chromflock_init();
    }

    if(!strcmp(command, "version"))
    {
        return show_version();
    }

    if(!strcmp(command, "hic2cpm"))
    {
        char * cmd = malloc(strlen(argv[0]) + strlen(argv[1]) + 2);
        sprintf(cmd, "%s %s", argv[0], argv[1]);
        argv[0] = cmd;
        for(int kk = 1; kk<argc; kk++)
        {
            argv[kk] = argv[kk+1];
        }
        int ret = cc2cpm(argc-1, argv);
        free(cmd);
        return ret;
    }

    if(!strcmp(command, "string2any"))
    {
        char * cmd = malloc(strlen(argv[0]) + strlen(argv[1]) + 2);
        sprintf(cmd, "%s %s", argv[0], argv[1]);
        argv[0] = cmd;
        for(int kk = 1; kk<argc; kk++)
        {
            argv[kk] = argv[kk+1];
        }
        int ret = string2any(argc-1, argv);
        free(cmd);
        return ret;
    }

    if(!strcmp(command, "any2string"))
    {
        char * cmd = malloc(strlen(argv[0]) + strlen(argv[1]) + 2);
        sprintf(cmd, "%s %s", argv[0], argv[1]);
        argv[0] = cmd;
        for(int kk = 1; kk<argc; kk++)
        {
            argv[kk] = argv[kk+1];
        }
        int ret = any2string(argc-1, argv);
        free(cmd);
        return ret;
    }

    if(!strcmp(command, "sprite2cpm"))
    {
        char * cmd = malloc(strlen(argv[0]) + strlen(argv[1]) + 2);
        sprintf(cmd, "%s %s", argv[0], argv[1]);
        argv[0] = cmd;
        for(int kk = 1; kk<argc; kk++)
        {
            argv[kk] = argv[kk+1];
        }
        int ret = sprite2cmap(argc-1, argv);
        free(cmd);
        return ret;
    }

    usage();
    return EXIT_FAILURE;
}
