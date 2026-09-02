#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include "collide3.h"

typedef uint32_t u32;
typedef uint64_t u64;

#ifdef COLLIDE3_f64
typedef double fxx;
#define fxx(X) X##_f64
#else
typedef float fxx;
#define fxx(X) X##_f32
#endif

// hash table entry
typedef struct {
    fxx X[3]; // coordinate
    u32 idx; // original index
} entry;

static fxx
eudist2(const fxx * A, const fxx * B)
{
    // Euclidean distance between two 3D-vectors to the power of two
    // i.e., ||A-B||^2
    return pow(A[0]-B[0], 2) + pow(A[1]-B[1], 2) + pow(A[2]-B[2], 2);
}

static double
timespec_diff(struct timespec* end, struct timespec * start)
{
    double elapsed = (end->tv_sec - start->tv_sec);
    elapsed += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}

// Hash function for a single coordinate
static u32
hash_coord(const int nDiv, fxx x)
{
    x < -1.0 ? x = -1.0 : 0 ;
    x > 1.0-1e-6 ? x = 1.0-1e-6 : 0 ;
    assert( ((x+1.0)/2.0 * nDiv) < nDiv );
    return (u32) ((x+1.0)/2.0 * nDiv) ;
}

// Hash function for a 3D point
static u32
hash(const u32 nDiv, const fxx * X)
{
    u32 h = hash_coord(nDiv, X[0]) +
        nDiv*hash_coord(nDiv, X[1]) +
        nDiv*nDiv*hash_coord(nDiv, X[2]);
    return h;
}

// The default dummy callback used when no callback function is
// passed.
static void
dummy_cb(__attribute__((unused)) u32 u,
              __attribute__((unused)) u32 v,
              __attribute__((unused)) double d2,
              __attribute__((unused)) void * data)
{
    return;
}

int
fxx(collide3)(const fxx * D, const u32 N, const fxx d,
             collide3_info * info,
             collide3_cb cb, void * cb_data)
{
    struct timespec t0, t1, t2, t3;
    if(info) {
        clock_gettime(CLOCK_REALTIME, &t0);
    }

    if(cb == NULL) {
        cb = dummy_cb;
        cb_data = NULL;
    }

    const fxx d2 = powf(d,2); // store squared distance for future use

    // Make a decision about the number of buckets
    // per dimensions. Since we only care about the
    // a cubic domain, it will be the same for all dimensions
    //
    // The current heuristic was found using tests where
    // N=1,000,000, r = 2.0/cbrt(n), some results:
    // cbrt(N/4) -> 1,315.89 s
    // cbrt(N/5) -> 1,208.33 s
    // cbrt(N/6) -> 1,322.54 s
    // cbrt(N/7) -> 1,507.23 s
    // cbct(N/8) -> 1,810.87 s
    int nDiv = cbrt(N/5);
    nDiv < 2 ? nDiv = 2 : 0;
    u32 n_buckets = nDiv*nDiv*nDiv;

    //
    // Count sort to create the hash table, HT.
    //

    // We will reference this array under several aliases and offsets,
    // use paper and pen to figure out!
    u32 * bucket_list = calloc(n_buckets+3, sizeof(u32));
    if(bucket_list == NULL) { return EXIT_FAILURE; }
    if(info) {
        info->mem_alloc = (n_buckets+3)*sizeof(u32);
    }
    u32 * bucket_size = bucket_list + 2;

    // Count how many elements that will fall into each bucket.
    for(u32 kk = 0; kk<N; kk++) {
        bucket_size[hash(nDiv, D+3*kk)]++;
    }

    // Integrate the list -- find the start position of each bucket
    u32 * bucket_writepos = bucket_list + 1;
    for(u32 kk = 1; kk<=n_buckets; kk++) {
        bucket_writepos[kk] = bucket_size[kk-1]+bucket_writepos[kk-1];
    }

    // Create the actual hash table and insert copies
    // of all points tagged with their original indices
    // i.e. sort the points according to their bucket
    //
    // Points are copied to avoid cache misses in the later scanning
    // phase.
    //
    // Accumulates the the write positions
    // so that in the end writepos[kk] is startpos[kk+1]
    entry * HT = malloc(N*sizeof(entry));
    if(info) {
        info->mem_alloc += N*sizeof(entry);
    }
    if(HT == NULL) {
        free(bucket_list);
        return EXIT_FAILURE;
    }

    for(u32 kk = 0; kk<N; kk++) {
        u32 h = hash(nDiv, D+3*kk);
        u32 ht_pos = bucket_writepos[h]; // hash table position
        bucket_writepos[h]++;
        HT[ht_pos].idx = kk;
        HT[ht_pos].X[0] = D[3*kk];
        HT[ht_pos].X[1] = D[3*kk+1];
        HT[ht_pos].X[2] = D[3*kk+2];
    }

    u32 * bucket_start = bucket_list;

    if(info) {
        clock_gettime(CLOCK_REALTIME, &t1);
        info->t_create_ms = 1000.0 * timespec_diff(&t1, &t0);
        clock_gettime(CLOCK_REALTIME, &t2);
    }

    // The loop is over the elements of the hash table,
    // i.e. not over the order that the elements were given
    // under the idea that this might reduce cache misses
    u64 counter = 0;
    for(u32 kk = 0; kk<N; kk++) {
        // Figure out which bins might contain a hit
        fxx deps = d;
        u32 ha_min = hash_coord(nDiv, HT[kk].X[0] - deps);
        u32 ha_max = hash_coord(nDiv, HT[kk].X[0] + deps);
        u32 hb_min = hash_coord(nDiv, HT[kk].X[1] - deps);
        u32 hb_max = hash_coord(nDiv, HT[kk].X[1] + deps);
        u32 hc_min = hash_coord(nDiv, HT[kk].X[2]); // will never scan backward in the last dimension
        u32 hc_max = hash_coord(nDiv, HT[kk].X[2] + deps);

        // Loop over the "worst" i.e. most strided dimension first
        for(u32 cc = hc_min; cc <= hc_max; cc++) {
            for(u32 bb = hb_min; bb <= hb_max; bb++) {

                // last bucket to visit
                u32 hash1 =
                    ha_max +
                    bb*nDiv +
                    cc*pow(nDiv,2);

                if(bucket_start[hash1+1]-1 < kk) continue;

                // first bucket to visit
                u32 hash0 =
                    ha_min +
                    bb*nDiv +
                    cc*pow(nDiv,2);

                // index of first element to compare with
                u32 ht_start = bucket_start[hash0];
                // we only compare to later elements in the table
                ht_start <= kk ? ht_start = kk+1 : 0;
                // index of last element to compare with
                u32 ht_end = bucket_start[hash1+1];
                for(u32 pp = ht_start; pp < ht_end; pp++) {
                    fxx pd2 = eudist2(HT[pp].X, HT[kk].X);
                    if(pd2 < d2) { // squared distances
                        cb(HT[pp].idx, HT[kk].idx, pd2, cb_data);
                        counter++;
                    }
                }
            }
        }
    }
    free(bucket_list);
    free(HT);
    if(info) {
        clock_gettime(CLOCK_REALTIME, &t3);
        info->t_scan_ms = 1000.0 * timespec_diff(&t3, &t2);
        info->n_collisions = counter;
    }
    return EXIT_SUCCESS; // == success
}
