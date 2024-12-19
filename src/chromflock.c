/**
 * @file chromflock.c
 * @author Erik Wernersson
 * @date 2020-2023
 */

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
#include "contact_pairs_io.h"

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
           "convert raw data dumps to text\n");
    printf("string2any\n\t"
           "write human readable to raw\n");
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
        return cc2cpm(argc-1, argv+1);
    }

    if(!strcmp(command, "string2any"))
    {
        return string2any(argc-1, argv+1);
    }

    if(!strcmp(command, "any2string"))
    {
        return any2string(argc-1, argv+1);
    }

    if(!strcmp(command, "sprite2cpm"))
    {
        return sprite2cmap(argc-1, argv+1);
    }

    if(!strcmp(command, "unittests"))
    {
        return contact_pairs_io_ut(argc-1, argv+1);
    }


    fprintf(stderr, "%s is an unknown command to me\n", command);

    usage();
    return EXIT_FAILURE;
}
