#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <sys/types.h>

// Enable to see some timings
// #define  DDICT_TIMINGS

#ifdef DDICT_TIMINGS
#include <time.h>
#endif

#include "ddict.h"

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

// Settings
// The storage for entries will be half that size
#define DDICT_INDEX_INITIAL_SIZE 32
// The index will be at least this factor
// larger than the number of inserted elements
#define DDICT_INDEX_FACTOR 3
// To rebuild the index is quite expensive. A large grow rate
// will in general be faster but require more memory.
#define DDICT_INDEX_GROWTH_RATE 4
// factor to grow entry storage with when full
#define DDICT_ENTRIES_GROWTH_RATE 2
// factor to grow key storage with when full
#define DDICT_KEY_STORAGE_GROWTH_RATE 2

#ifdef DDICT_TIMINGS
static double
timespec_diff(struct timespec* end, struct timespec * start)
{
    double elapsed = (end->tv_sec - start->tv_sec);
    elapsed += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}
#endif


static char *
ddict_store_key(ddict * dict, const char *str)
{
    size_t nb = strlen(str) + 1;
    if(nb + dict->key_storage_pos >= dict->key_storage_size)
    {
        dict->key_storage_size *= DDICT_KEY_STORAGE_GROWTH_RATE;
        #ifdef DDICT_TIMINGS
        struct timespec t0, t1;
        clock_gettime(CLOCK_REALTIME, &t0);
        #endif
        u8 * new = realloc(dict->key_storage, dict->key_storage_size);
        if(new != dict->key_storage) {
            ssize_t delta = new - dict->key_storage;
            for(u64 kk = 1; kk <= dict->n_entries; kk++)
            {
                dict->entries[kk].key += delta;
            }
            dict->key_storage = new;
        }
#ifdef DDICT_TIMINGS
        clock_gettime(CLOCK_REALTIME, &t1);
        printf("Increase key storage capacity took %f s\n", timespec_diff(&t1, &t0));
#endif
    }
    memcpy(dict->key_storage + dict->key_storage_pos, str, nb);
    char * retpos = (char*) (dict->key_storage + dict->key_storage_pos);
    dict->key_storage_pos += nb;
    #if 0
    for(u64 kk = 0; kk < dict->key_storage_pos; kk++) {
        char c = dict->key_storage[kk];
        if(c == '\0') {
            printf("|");
        } else {
            printf("%c", c);
        }
    }
    printf("\n");
    #endif
    return retpos;
}

#if 0
static u64
_ddict_wordhash(const char * word)
{
    return (u64) word;
}
#endif

static u64
_ddict_wordhash(const char * word)
{
    // your standard Polynomial rolling hash function
    u8 * bytes = (u8*) word;
    u64 h = 0;
    while(*bytes != '\0') {
        h = h*31 + *bytes++;
    }
    return h;
}

static inline u64
_ddict_entry_id_from_index(const ddict * dict, const u64 idx)
{
    switch(dict->index_type)
    {
    case U8:
        return (u64) dict->index.U8[idx];
    case U16:
        return (u64) dict->index.U16[idx];
    case U32:
        return (u64) dict->index.U32[idx];
    case U64:
        return (u64) dict->index.U32[idx];
    }
    __builtin_unreachable();
    assert(0);
    return -1;
}

static inline u64
#ifdef DDICT_STATS
_ddict_get_free_index(ddict * dict, u64 idx)
#else
    _ddict_get_free_index(const ddict * dict, u64 idx)
#endif
{
    switch(dict->index_type)
    {
    case U8:
        while(dict->index.U8[idx] != 0) {
            idx++;
            if(idx == dict->n_indices) {
                idx = 0;
            }
#ifdef DDICT_STATS
            dict->n_collision++;
#endif
        }
        break;
    case U16:
        while(dict->index.U16[idx] != 0) {
            idx++;
            if(idx == dict->n_indices) {
                idx = 0;
            }
#ifdef DDICT_STATS
            dict->n_collision++;
#endif
        }
        break;
    case U32:
        while(dict->index.U32[idx] != 0) {
            idx++;
            if(idx == dict->n_indices) {
                idx = 0;
            }
#ifdef DDICT_STATS
            dict->n_collision++;
#endif
        }
        break;
    case U64:
        while(dict->index.U32[idx] != 0) {
            idx++;
            if(idx == dict->n_indices) {
                idx = 0;
            }
#ifdef DDICT_STATS
            dict->n_collision++;
#endif
            break;
        }
    }
    return idx;
}

static inline void
_ddict_set_index(ddict * dict, const u64 idx, const u64 value)
{
    assert(value != 0);
    switch(dict->index_type) {
    case U8:
        dict->index.U8[idx] = value;
        return;
    case U16:
        dict->index.U16[idx] = value;
        return;
    case U32:
        dict->index.U32[idx] = value;
        return;
    case U64:
        dict->index.U32[idx] = value;
        return;
    }
}

// Create the indices array, selecting the appropriate
// element size depending on the number of indices
static int
_ddict_gen_indices(ddict * dict) {
    // select data type
    dict->index_type = U8;
    size_t esize = 1;
    if(dict->n_indices > 256-2){
        dict->index_type = U16;
        esize = 2;
    }
    if(dict->n_indices > 65025-2){
        dict->index_type = U32;
        esize = 4;
    }
    if(dict->n_indices > 4294967296-2){
        dict->index_type = U64;
        esize = 8;
    }

    dict->index.U8 = calloc(dict->n_indices, esize);
    if(dict->index.U8 == NULL) {
        return 1;
    }

    return 0;
}

static ddict *
ddict_new_with_size(const int manage_keys,
                    const u64 n_indices,
                    const int alloc_entries)
{

    if(n_indices < 2) {
        return NULL;
    }

    ddict * dict = calloc(1, sizeof(ddict));
    if(dict == NULL) {
        return NULL;
    }
    dict->manage_keys = manage_keys;
    dict->n_indices = n_indices;
    if(_ddict_gen_indices(dict)) {
        free(dict);
        return NULL;
    }

    if(!alloc_entries) {
        return dict;
    }

    if(manage_keys)
    {
        dict->key_storage_size = 20*n_indices;
        dict->key_storage = malloc(dict->key_storage_size);
        if(dict->key_storage == NULL) {
            free(dict);
            return NULL;
        }
    }

    dict->n_entries_alloc = n_indices / 2;
    dict->entries = malloc((dict->n_entries_alloc+1)*sizeof(ddict_entry));
    if(dict->entries == NULL) {
        free(dict->index.U8);
        free(dict->key_storage);
        free(dict);
        return NULL;
    }

    return dict;
}

ddict *
ddict_new(const int manage_keys) {
    return ddict_new_with_size(manage_keys, DDICT_INDEX_INITIAL_SIZE, 1);
}

void
ddict_free(ddict * dict) {
    if(dict == NULL) {
        return;
    }
#ifdef DDICT_STATS
    printf("ddict_free (#DDICT_STATS defined)\n"
           "            dict->n_collision = %lu\n"
           "            dict->n_indices = %lu\n",
           dict->n_collision, dict->n_indices);
#endif
    free(dict->key_storage);
    free(dict->entries);
    free(dict->index.U8);
    free(dict);
    return;
}


static ddict_entry *
ddict_get_with_hash(const ddict * dict,
                    const char * key, const u64 hash)
{
    u64 idx = hash % dict->n_indices;

    // Until found or an empty slot appears
    while(1)
    {
        if(idx == dict->n_indices) {
            idx = 0; // wrap around
        }

        // i64 eid = dict->indices[idx]; // ddict_entry index
        u64 eid = _ddict_entry_id_from_index(dict, idx);
        if(eid == 0) {
            return NULL;
        }

#ifndef DDICT_DROP_HASH
        if(dict->entries[eid].hash != hash)
        {
            idx++;
            continue;
        }
#endif

        if(strcmp(key, dict->entries[eid].key) == 0)
        {
            return &(dict->entries[eid]);
        }
        idx++;
    }
    assert(0);
    __builtin_unreachable();
    return NULL;
}

ddict_entry *
ddict_get(const ddict * dict,
          const char * key)
{
    const u64 hash = _ddict_wordhash(key);
    return ddict_get_with_hash(dict, key, hash);
}


// Make room for more entries, This can have a separate growth rate
// compared to the index as it does not matter if it is almost full.
static int
ddict_grow_entries(ddict * dict)
{
    dict->n_entries_alloc = DDICT_ENTRIES_GROWTH_RATE * dict->n_entries_alloc + 1;
    dict->entries = realloc(dict->entries,
                            dict->n_entries_alloc*sizeof(ddict_entry));
    if(dict->entries == NULL)
    {
        return -1;
    }
    return 0;
}


static int
ddict_grow_indices(ddict * dict)
{
    //printf("!!! Growing indices (%lu)\n", dict->n_indices);
    // Procedure
    // - create a new dict (dict2) that shares the entries from the original dict
    // - in dict2, index all entries present in the original dict
    // - transfer the index in dict2 to the original dict

#ifdef DDICT_TIMINGS
    struct timespec t0, t1;
    clock_gettime(CLOCK_REALTIME, &t0);
#endif
    assert(dict->n_indices > 0);
    u64 n_indices = dict->n_indices;
    u64 n_indices2 = n_indices*DDICT_INDEX_GROWTH_RATE;
    free(dict->index.U8);
    ddict * dict2 = ddict_new_with_size(dict->manage_keys, n_indices2, 0);
    if(dict2 == NULL) { return -1; }
    dict2->entries = dict->entries;
    dict2->n_entries_alloc = dict->n_entries_alloc;

    for(u64 kk = 1; kk <= dict->n_entries; kk++)
    {
        ddict_entry e = dict->entries[kk];
#ifdef DDICT_DROP_HASH
        u64 idx = _ddict_wordhash(e.key) % n_indices2;
#else
        u64 idx = e.hash % n_indices2;
#endif
        // Find a free index slot
        idx = _ddict_get_free_index(dict2, idx);
        // Make it point to the already present entry
        _ddict_set_index(dict2, idx, kk);
    }

    // Note: No need to switch since all point to the same memory location
    dict->index.U32 = dict2->index.U32;
    dict->index_type = dict2->index_type;

#ifdef DDICT_STATS
    dict->n_collision = dict2->n_collision;
#endif
    dict->n_indices = dict2->n_indices;
    // We steal the indices from dict2 and it never had an entry
    // so we can just free the struct
    free(dict2);
#ifdef DDICT_TIMINGS
    clock_gettime(CLOCK_REALTIME, &t1);
    printf("grow indices took %f s\n", timespec_diff(&t1, &t0));
#endif

    return 0;
}


int
ddict_add(ddict * dict,
          char * word, void * value)
{
    const u64 hash = _ddict_wordhash(word);

    // See if in dictionary
    if(ddict_get_with_hash(dict, word, hash) != NULL) {
        return 1;
    }

    // Increase the entry capacity if needed
    if(dict->n_entries+1 == dict->n_entries_alloc) {
        if(ddict_grow_entries(dict)){
            return -1;
        }
    }

    // Increase the index size of needed
    if(DDICT_INDEX_FACTOR*dict->n_entries > dict->n_indices) {
        if(ddict_grow_indices(dict)) {
            return -1;
        }
    }

    // Append the new entry to entries
#ifndef DDICT_DROP_HASH
    dict->entries[dict->n_entries+1].hash = hash;
#endif
    if(dict->manage_keys) {
        dict->entries[dict->n_entries+1].key = ddict_store_key(dict, word);
    } else {
        dict->entries[dict->n_entries+1].key = word;
    }

    dict->entries[dict->n_entries+1].value = value;

    // Figure out where we can insert a reference
    u64 idx = hash % dict->n_indices;
    idx = _ddict_get_free_index(dict, idx);
    _ddict_set_index(dict, idx, dict->n_entries+1);
    //printf("idx = %lu\n", idx);
    dict->n_entries++;
    return 0;
}

u64
ddict_size(const ddict * dict)
{
    return dict->n_entries;
}

int ddict_update_entry(ddict * dict, const char * key, void * value)
{
    ddict_entry * e;
    if( (e = ddict_get(dict, key)) == NULL)
    {
        return 1;
    }
    e->value = value;
    return 0;
}
