#include "chromflock_init.h"


/* Where is the binary located? Returns "./" if it can't figure out. */
char * mylocation(void)
{
    size_t bufsize = 1024;
    char * buf = malloc(bufsize);
    if(readlink("/proc/self/exe", buf, bufsize) == -1)
    {
        sprintf(buf, "./");
    }
    return buf;
}

int isfile(const char * filename)
{
    struct stat statbuf;
    if(stat(filename, &statbuf))
    {
        return 0;
    }

    /* see man inode(7) */
    if( !(S_ISREG(statbuf.st_mode) ))
    {
        return 0;
    }
    return 1;
}

/* Convert
   /home/user/dir/chromflock/bin/chromflock to
   /home/user/dir/chromflock/util/chromflock_gen
   Returns "/chromflock_gen" if that isn't possible
*/
char * get_src_dir_chromflock_gen(void)
{
    char * filename = malloc(2048);
    char * myloc = mylocation();
    sprintf(filename, "%s", "/chromflock_gen");

    /* If at least three '/' then replace the second last with '\0'' */
    int nslash = 0;
    for(size_t kk = 0; kk < strlen(myloc); kk++)
    {
        if( myloc[kk] == '/')
        {
            nslash++;
        }
    }

    if(nslash > 2)
    {
        int n = 0;
        for(size_t kk = 0; kk < strlen(myloc); kk++)
        {
            if( myloc[kk] == '/')
            {
                n++;
                if(n == nslash - 1)
                {
                    myloc[kk] = '\0';
                }
            }
        }

        sprintf(filename, "%s/util/chromflock_gen", myloc);
    }
    free(myloc);
    return filename;
}

int chromflock_init(void)
{
    /* find the chromflock_gen script
     *  and copy it to the current folder
     */


    char * outname = strdup("chromflock_gen");
    char * filename = NULL;

    if(isfile(outname))
    {
        fprintf(stderr,
                "A file called '%s' already exists in this folder "
                "please remove it before using this command\n",
                outname);
        goto fail;
    }

    filename = get_src_dir_chromflock_gen();
    if(isfile(filename))
    {
        goto copy_file;
    } else {
        printf("%s does not exist\n", filename);
    }

    free(filename);
    filename = strdup("/usr/local/share/chromflock/chromflock_gen");
    if(isfile(filename))
    {
        goto copy_file;
    } else {
        printf("%s does not exist\n", filename);
    }

    fprintf(stderr,
            "Unable to find the file chromflock_gen! This indicates that the "
            "installation is incomplete. Please make sure that chromflock is "
            "properly installed.\n");
    goto fail;

copy_file:

    printf("Found source script at: %s\n", filename);
    if(oscp(filename, outname) <= 0)
    {
        fprintf(stderr,
                "Failed to copy %s to %s. Please check that you have read "
                "permission for the former and that you have write "
                "permissions in the current folder\n",
                filename, outname);
        goto fail;
    }

    printf("Created the file '%s' in this folder.\n"
           "Please edit it for your needs.\n"
           "When done, run:\n"
           "  chmod +x %s\n"
           "  ./%s\n",
           outname, outname, outname);

    free(filename);
    free(outname);
    return EXIT_SUCCESS;

fail:
    free(filename);
    free(outname);
    return EXIT_FAILURE;
}
