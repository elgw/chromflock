#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdlib.h>

// A collision detector for points in [-1, 1]^3
//
// Implemented using spatial hashing and count sort.
//
// Version 1.0.0
// Erik Wernersson 20260901

    typedef void (*collide3_cb)(
        uint32_t u, // the index of one point
        uint32_t v, // the index of the other point
        double d2, // the squared Euclidean distance between the points
        void * cb_data); // user pointer

    typedef struct {
        double t_create_ms;
        double t_scan_ms;
        uint32_t mem_alloc; // number of allocated bytes
        uint64_t n_collisions;
    } collide3_info;


// Detect collisions among the 3D points, i.e., where ||u-v|| < radius
//
// Arguments:
// - points: [[x, y, z], [x, y, z], ...].
//
// - cb: a function to be called upon a collision detection. Will only
//   be called once per pair, i.e. if the points u, v are in contact,
//   the function will be called for for either (u, v) or (v, u). Can
//   be NULL.
//
// - info The fields of this struct will be populated. Can be NULL.
//
// Returns EXIT_SUCCESS unless a memory allocation error occured.
//
// Important: The methods does not check that the points are within
// the bounds. Points outside of the bounds will cause INEFFICIENT
// queries.


    int
    collide3_f32(const float * points,
                 uint32_t n_point,
                 float point_radius,
                 collide3_info * info,
                 collide3_cb cb,
                 void * cb_data);

    // Only available if COLLIDE3_f64 is defined
    // at build
    int
    collide3_f64(const double * points,
                 uint32_t n_point,
                 double point_radius,
                 collide3_info * info,
                 collide3_cb cb,
                 void * cb_data);


#ifdef __cplusplus
}
#endif
