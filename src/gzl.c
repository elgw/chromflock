#include <stdio.h>
#include <string.h>

#include "gzl.h"

typedef int64_t i64;

struct gzl__state {
    gzFile zf; /* zlib handle to open file */
    uint8_t * buf; /* buffer for storing read data */
    int64_t buf_size; /* size of buffer, has to be at least the line size */
    int64_t buf_used; /* number of bytes used */
    int64_t next_pos; /* where the next line starts in the buffer */
    const char * error_string;
    int error_code;
};

gzl_state * gzl_open(const char * filename, i64 max_line_size)
{
    gzFile zf = gzopen(filename, "rb");
    if(zf == NULL)
    {
        return NULL;
    }
    gzl_state * gzl = calloc(1, sizeof(gzl_state));
    gzl->buf_size = max_line_size;
    // The buffer will contain at last one line
    // One byte reserved for writing \0
    gzl->buf = calloc(1, gzl->buf_size + 1);
    gzl->buf_used = 0;
    gzl->zf = zf;
    gzl->error_string = "0 is not an error";
    return gzl;
}

void gzl_destroy(gzl_state * gzl)
{
    free(gzl->buf);
    gzclose(gzl->zf);
    free(gzl);
    return;
}

const char * gzl_get_line(gzl_state * gzl, int * error)
{
    if(gzl->error_code != 0)
    {
        return NULL;
    }
    /* Fast path: Find another line in the buffer and return it. */
    for(i64 kk = gzl->next_pos; kk < gzl->buf_used; kk++)
    {
        if( (gzl->buf[kk] == '\n') | (gzl->buf[kk] == '\0'))
        {
            gzl->buf[kk] = '\0';
            const char * str = (const char *) gzl->buf + gzl->next_pos;
            gzl->next_pos = kk+1;
            return str;
        }
    }

    /* Read more data
       - Cut out last line: i.e. shift remaining data
       - Fill the buffer
    */
    if(gzl->buf_used >= gzl->next_pos)
    {
        i64 nmove = gzl->buf_used - gzl->next_pos;
        memmove(gzl->buf,
                gzl->buf + gzl->next_pos,
                nmove);

        gzl->buf_used = nmove;
        gzl->next_pos = 0;
    }


    /* Read more data */
    if(gzl->buf_size - gzl->buf_used == 0)
    {
        gzl->error_code = 1;
        gzl->error_string  = "Found a line longer than the max_line_size";
        *error = 1;
        return NULL;
    }

    if(gzl->buf_size - gzl->buf_used < 0)
    {
        gzl->error_code = 2;
        gzl->error_string =
            "internal library error, something unexpected happened, "
            "please file a bug report!";
        *error = 2;
        return NULL;
    }
    int status = gzread(gzl->zf,
                        gzl->buf + gzl->buf_used,
                        gzl->buf_size - gzl->buf_used);
    if(status <= 0)
    {
        goto end_or_err;
    }
    gzl->buf_used += status;

    uint8_t * p = gzl->buf;
    i64 pos = 0;
    while( (pos < gzl->buf_size) && (p[pos] != '\n'))
    {
        pos++;
    }
    gzl->next_pos = pos+1;
    gzl->buf[pos] = '\0';
    return (const char *) gzl->buf;

end_or_err:
    ;
    // libz uses negative values for errors
    int zerrno;
    gzl->error_string = gzerror(gzl->zf, &zerrno);
    if(zerrno == 0)
    {  // End of file
        return NULL;
    }
    *error = zerrno;
    return NULL;
}

const char * gzl_get_error_string(gzl_state * gzl)
{
    return gzl->error_string;
}
