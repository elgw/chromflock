#include "chromflock_init.h"


/* Where is the binary located? Returns "./" if it can't figure out. */
char * mylocation(void)
{
    size_t bufsize = 1024;
    char * buf = malloc(bufsize);
    assert(buf != NULL);

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
   /home/user/dir/chromflock/
   Returns "/" if that isn't possible
*/
char * get_src_dir(void)
{
    char * myloc = mylocation();

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
    }
    return myloc;
}

typedef struct {
    char * filename;
    char * relpath;
    mode_t mode;
} cf_file;


/* Copy the file, specified by name to this folder.
 * If chromflock is run from the source dir in relpath, prefer files in the
 * repository. Else look in the default install dir. */
int place_here(const cf_file * file)
{

    char * outname = file->filename;

    // Check if file already exist
    if(isfile(outname))
    {
        fprintf(stderr,
                "ERROR: A file called '%s' already exists in this folder "
                "please remove it before using this command\n",
                outname);
        return EXIT_FAILURE;
        goto fail;
    }

    // Copy from source directory
    char * srcdir = get_src_dir();
    char * filename = malloc( strlen(srcdir)
                              +strlen(file->relpath)
                              +strlen(file->filename)+16 );
    assert(filename != NULL);
    sprintf(filename, "%s%s%s", srcdir, file->relpath, file->filename);
    free(srcdir);
    if(isfile(filename))
    {
        goto copy_file;
    }

    sprintf(filename, "/usr/local/share/chromflock/%s", file->filename);

    if(isfile(filename))
    {
        goto copy_file;
    }

    fprintf(stderr,
            "ERROR: Unable to find the file %s! This indicates that the "
            "installation is incomplete. Please make sure that chromflock is "
            "properly installed.\n", file->filename);
    goto fail;

copy_file:

    printf("Found source script at: %s\n", filename);
    if(oscp(filename, outname) <= 0)
    {
        fprintf(stderr,
                "ERROR: Failed to copy %s to %s. Please check that you have read "
                "permission for the former and that you have write "
                "permissions in the current folder\n",
                filename, outname);
        goto fail;
    }

    if(chmod(outname, file->mode))
    {
        printf("Failed to change mode for %s\n", file->filename);
    }

    free(filename);

    return EXIT_SUCCESS;

fail:
    free(filename);
    return EXIT_FAILURE;
}


int chromflock_init(void)
{

    /* Copy the following files to the current directory
     * */

    cf_file files[] = {
        {"chromflock_gen", "/util/", S_IRUSR | S_IWUSR | S_IXUSR },
        {"mflock.lua", "/src/", S_IRUSR | S_IWUSR},
        {NULL, NULL, 0}
    };

    int ok = 1;
    cf_file * file = &files[0];
    while( file->filename != NULL )
    {
        if(place_here(file))
        {
            ok = 0;
        }
        file++;
    }

    if(ok)
    {
        printf("Created chromflock_gen and mflock.lua in this folder.\n"
               "Edit them according to your needs.\n"
               "When done, run:\n"
               "  ./chromflock_gen\n");

    } else {
        printf("Failed to initialize this directory for chromflock. "
               "Please see previous error messages\n");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
