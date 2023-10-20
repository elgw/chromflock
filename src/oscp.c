#include "oscp.h"

int oscp(const char * source, const char * destination)
{
    int input, output;

    if ((input = open(source, O_RDONLY)) == -1)
    {
        return EXIT_FAILURE;
    }

    if ((output = creat(destination, 0660)) == -1)
    {
        close(input);
        return EXIT_FAILURE;
    }

#if defined(__APPLE__) || defined(__FreeBSD__)
    //fcopyfile works on FreeBSD and OS X 10.5+
    int result = fcopyfile(input, output, 0, COPYFILE_ALL);
    if(result < 0)
    {
        fprintf(stderr, "Failed to copy %s to %s\n", source, destination);
        close(input);
        close(output);
        return EXIT_FAILURE;
    }
#else
    //sendfile will work with non-socket output (i.e. regular file) on Linux 2.6.33+
    off_t bytesCopied = 0;
    struct stat fileinfo = {0};
    fstat(input, &fileinfo);
    ssize_t result = sendfile(output, input, &bytesCopied, fileinfo.st_size);
    if(result == -1)
    {
        fprintf(stderr, "Failed to copy %s to %s\n", source, destination);
        close(input);
        close(output);
        return EXIT_FAILURE;
    }
#endif

    close(input);
    close(output);

    return result;
}
