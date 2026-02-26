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
#include "con2mflock.h"
#include "chromflock_usage.h"

static int show_version(void)
{
    printf("chromflock v.%s\n", cf_version);
    return EXIT_SUCCESS;
}

static int usage(void)
{
    printf("%s", src_txt_chromflock_usage_txt);
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

    if(!strcmp(command, "con2mflock"))
    {
        return con2mflock(argc-1, argv+1);
    }

    fprintf(stderr, "%s is an unknown command to me\n", command);

    usage();
    return EXIT_FAILURE;
}
