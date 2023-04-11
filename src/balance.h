#ifndef __balance_h__
#define __balance_h__

/* Matrix balancing of A so that sum(M,1) = sum(M,2)' = 1
 * Columns/rows that sum to 0 are ignored.
 * Assumes that M^T=M
 * Assumes that all values of M are normal/finite.
 * Returns the max absolute difference between any
 * column sum and 1, while ignoring columns that are all zeros.
 * Returns -1 in the case that M is 0
 */

#include "stdint.h"

double balance(double * M, size_t N);

#endif
