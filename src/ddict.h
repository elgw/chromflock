#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdint.h>

// A dense key storage dictionary data structure, similar to what is
// found in Python but less sophisticated
//
// - Grows dynamically when needed
//
// - Using Open Addressing with Linear Probing
//
// - Can either own all the keys, i.e. makes copies of them and make
//   sure to free them at the end, or have the caller own the
//   keys. This is controlled by the argument to ddict_new()
//
// - Only string keys, will be scanned until a '\0' is found
//
// - If the hash value is cheap to compute, there is no need to save
//   them, disable hash storage by defining DDICT_DROP_HASH below.
//
// Erik Wernersson, 2026-08-14

    #define DDICT_VERSION_MAJOR 1
    #define DDICT_VERSION_MINOR 0
    #define DDICT_VERSION_PATCH 2

// If DDICT_DROP_HASH is defined, the hash values will not be
// calculated when needed, i.e. not stored.
// #define DDICT_DROP_HASH

// Records the total number of steps the insertions happened
// away from the ideal locations. Useful to check when switching
// hash function
// #define DDICT_STATS

    typedef struct {
#ifndef DDDICT_DROP_HASH
        uint64_t hash;
#endif
        char * key;
        void * value;
    } ddict_entry;

// Memory layout
//
//                      First entry is always empty
//                      |
// ddict->entries = [ [0,0,0], [ 9387453, "john", "doe"], [36509, "homer", "simpson"], ... ]
// ddict->indices = [ 0 , 0, 2, 0, 1, 0, 0, 0]
//                    |      |     |
//                    |      |     location of john
//                    |      location of homer
//                    A zero means that the there is no entry associated
//
// The HASH value, h, points to a location in e=ddict->indices[h]
// If e is 0, there is nothing stored with that hash value. If non-zero, check
// ddict->entries[e], e_k=ddict->entries[ddict->indices[h+k]], etc until e_k is
// is 0 or a match is found.

    typedef enum {U8, U16, U32, U64} indices_bits;

    typedef union { // anonymous unions are not part of C99 :(
        uint8_t * U8;
        uint16_t * U16;
        uint32_t * U32;
        uint64_t * U64;
    } index_union;

    typedef struct {
        // Array of indices that points into the entries
        // or has the value 0 if there is nothing to be found
        indices_bits index_type;
        index_union index;

        // Size/number of elements of the index.
        uint64_t n_indices;

        // Storage for (hash, key, value) triplets. The first index is left unused.
        ddict_entry * entries;
        // Number of entries that are used
        uint64_t n_entries;
        // Number of entries that can be used (the actual allocation is
        // one more element).
        uint64_t n_entries_alloc;

        // If the dict allocate and store private copies
        // of the keys. The values are never owned the dict.
        int manage_keys;

        uint8_t * key_storage;
        uint64_t key_storage_size;
        uint64_t key_storage_pos; // where to write

#ifdef DDICT_STATS
        // Count the number of collisions during ddict_add
        uint64_t n_collision;
#endif

    } ddict;

// Create a new dictionary with the default size
    ddict * ddict_new(int manage_keys);

// Free the dictionary and all the key copies that it holds
    void ddict_free(ddict * dict);

// Return an ddict_entry with the given key or NULL if
// nothing is found
    ddict_entry * ddict_get(const ddict *, const char * key);

// Add an ddict_entry to the dictionary unless it already contains
// the given key. See also ddict_update_entry.
//
// Returns 0 on success.
    int
    ddict_add(ddict * dict,
              char * key, void * value);

// Return the number of entries in the dictionary
    uint64_t
    ddict_size(const ddict * dict);

// Update an existing key to point to a new value
    int
    ddict_update_entry(ddict * dict, const char * key, void * value);

#ifdef __cplusplus
}
#endif
