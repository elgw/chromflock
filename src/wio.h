// Compile with -lz

// Read an array stored as uint8_t. Use lz if the file ends with '.gz'
// returns NULL on failure
// WARNING: uses the last 4 bytes of 'gz' files to determine file size
// that fail if the file is  2^32 bytes or larger uncompressed.
void * wio_read(char * fileName, size_t * nel);

// Write an array stored as uint8_t. Use lz if the file ends with '.gz'
// returns 0 on sucess
int wio_write(char * fileName, size_t nel, void * data);
