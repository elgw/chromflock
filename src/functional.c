#include "functional.h"

static double norm3(const double * restrict X)
{
    double n = 0;
    for(size_t kk = 0; kk<3; kk++)
        n+=pow(X[kk], 2);
    return sqrt(n);
}

static void vec3_normalize(double * restrict X)
{
    double n = 0;
    for(size_t kk = 0; kk<3; kk++)
        n+=pow(X[kk], 2);
    X[0] /= n; X[1] /= n; X[2] /= n;
}


static double norm32(const double * restrict X)
{
    double n = 0;
    for(size_t kk = 0; kk<3; kk++)
        n+=pow(X[kk], 2);
    return n;
}

static double eudist3sq(const double * A, const double * B)
{
    /* SQUARED Euclidean distance between two 3D-vectors */
    return pow(A[0]-B[0], 2) + pow(A[1]-B[1], 2) + pow(A[2]-B[2],2);
}

static double eudist3(const double * A, const double * B)
{
    /* Euclidean distance between two 3D-vectors */
    return sqrt( pow(A[0]-B[0], 2) + pow(A[1]-B[1], 2) + pow(A[2]-B[2], 2));
}

static size_t hash_coord(const int nDiv, const double X)
{
    double v = (X+1)/2 * nDiv;

    #ifndef NDEBUG
    if(isfinite(v) != 1)
    {
        fprintf(stderr, "ERROR in %s, line %d\n", __FILE__, __LINE__);
        fprintf(stderr, "X=%f, v = %f\n", X, v);
        exit(EXIT_FAILURE);
    }
    #endif

    if(v<=0)
        return 0;
    if(v>=nDiv)
        return nDiv-1;
    return floor(v);
}


static size_t hash(const size_t nDiv, const double * restrict X)
{

    size_t h = hash_coord(nDiv, X[0]) +
        nDiv*hash_coord(nDiv, X[1]) +
        nDiv*nDiv*hash_coord(nDiv, X[2]);

    return h;
}


double errRepulsion(const double * restrict D,
                    const size_t N,
                    const double d)
{

    // D: coordinates of the dots [d0_x, d0_y, d0_z, d1_x, d1_y, ... dN-1_z]
    // N: Number of dots
    // d: capture distance between bead centers

    // NOTE: Will crash if D contains nan or inf values.

#ifndef NDEBUG
    // Verify that D was actually allocated
    for(size_t kk = 0; kk<3*N; kk++)
    {
        assert(isfinite(D[kk]));
    }
    fflush(stdout);
#endif


    // Mostly copied from volHASH/src/volhash3.c
    double d2 = pow(d,2);

    // 1. Counting -- Figure out how many elements per bucket
    int nDiv = 9; // Will create nDiv^3 buckets over [-1,1]^3

    if(N>10000)
    {
        nDiv = 15;
    }

    // i.e. bucket side length is 2/nDiv which can be related to d.
    // You can figure out optimal nDiv vs d and N by
    // running functinal_optimal_nDiv (while changing nDiv here)
    // For 3300 points: 9 or 10
    // For 33000 points: 15 about 30% faster than using 9

    size_t nH = pow(nDiv, 3);
    uint32_t * S = calloc(nH, sizeof(uint32_t));
    assert(S != NULL);

    for(size_t kk = 0; kk<N; kk++)
    {
        S[hash(nDiv, D+3*kk)]++;
    }

    /* Create boundaries for the data */
    uint32_t * B = malloc((nH+1)*sizeof(uint32_t));
    assert(B != NULL);
    uint32_t * C = malloc((nH+1)*sizeof(uint32_t));
    assert(C != NULL);

    B[0] = 0;
    C[0] = 0; // Start positions when writing
    for(size_t kk = 1; kk<=nH; kk++)
    {
        B[kk] = B[kk-1]+S[kk-1];
        C[kk] = B[kk];
    }

    if(0){
        for(size_t kk = 0; kk<nH; kk++)
        {
            printf("[%u -- %u], S[%zu] = %u\n", B[kk], B[kk+1], kk, S[kk]);
        }
    }

    /* Dots sorted according to their bucket and put into E */
    double * E = malloc(3*N*sizeof(double));
    assert(E!=NULL);

    /* Move data into new structure */
    for(size_t kk = 0; kk<N; kk++)
    {
        size_t h = hash(nDiv, D+3*kk);

        size_t writepos = 3*C[h];
        C[h]++;
        for(int idx = 0; idx<3; idx++)
        {
            assert(writepos<N*3);
            E[writepos++] = D[3*kk+idx];
        }
    }

    /*
     * now we can hash(X) a point X,
     * using B to see where the bucket is
     * in E
     * */

    // Get some job done!
    double err = 0;
    for(size_t kk = 0; kk<N; kk++)
        /* Find buckets that might contains neighbours */
    {
        double deps = d;
        size_t ha_min = hash_coord(nDiv, E[3*kk]-deps);
        size_t ha_max = hash_coord(nDiv, E[3*kk]+deps);

        size_t hb_min = hash_coord(nDiv, E[3*kk+1]-deps);
        size_t hb_max = hash_coord(nDiv, E[3*kk+1]+deps);

        size_t hc_min = hash_coord(nDiv, E[3*kk+2]-deps);
        size_t hc_max = hash_coord(nDiv, E[3*kk+2]+deps);

        for(size_t cc = hc_min; cc <= hc_max; cc++)
            for(size_t bb = hb_min; bb <= hb_max; bb++)
                for(size_t aa = ha_min; aa <= ha_max; aa++)
                {

                    // hash or index of the bucket to compare against
                    size_t hash =
                        aa +
                        bb*nDiv +
                        cc*pow(nDiv,2);

                    //          printf("%d, %d, %d, hash: %zu / %zu\n", aa, bb, cc, hash, nH); fflush(stdout);

                    // Compare against all elements in a specific bucket.
                    for(size_t pp = B[hash]; pp<B[hash+1]; pp++)
                    {
                        if(pp>kk) {
                            double dist2 = eudist3sq(E+3*pp, E+3*kk); //  squared distance
                            if(dist2<d2)
                            {
                                err += pow(sqrt(dist2) - d, 2);
                            }
                        }
                    }
                }
    }

    free(B);
    free(C);
    free(E);
    free(S);

    return err;
}


/** @brief Repulsion gradient (volumetric overlap)
 *
 * @param D The dots [3 X N]
 * @param G The gradient
 * @param N the number of dots
 * @param d ???
 * @param kVol The force magnitude
*/

static double
gradRepulsion(const double * restrict D,
                            double * restrict G,
                            const size_t N,
                            const double d,
                            const double kVol)
{

    assert(isfinite(d));

#ifndef NDEBUG
    // Verify that D was actually allocated
    for(size_t kk = 0; kk<3*N; kk++)
    {
        assert(isfinite(D[kk]));
    }
    fflush(stdout);
#endif

    // Mostly copied from volHASH/src/volhash3.c
    double d2 = pow(d,2);

    // 1. Figure out how many elements per bucket
    int nDiv = 9;
    if(N>10000)
    {
        nDiv = 15;
    }
    size_t nH = pow(nDiv, 3);
    uint32_t * S = calloc(nH, sizeof(uint32_t));
    assert(S!=NULL);

    for(size_t kk = 0; kk<N; kk++)
    {
        S[hash(nDiv, D+3*kk)]++;
    }

    /* Create boundaries for the data */
    uint32_t * B = malloc((nH+1)*sizeof(uint32_t));
    assert(B!=NULL);
    uint32_t * C = malloc((nH+1)*sizeof(uint32_t));
    assert(C!=NULL);

    B[0] = 0;
    C[0] = 0; // Start positions when writing
    for(size_t kk = 1; kk<=nH; kk++)
    {
        B[kk] = B[kk-1]+S[kk-1];
        C[kk] = B[kk];
    }

    if(0){
        for(size_t kk = 0; kk<nH; kk++)
        {
            printf("[%u -- %u], S[%zu] = %u\n", B[kk], B[kk+1], kk, S[kk]);
        }
    }

    // Dots sorted according to their bucket
    double * E = malloc(3*N*sizeof(double));
    assert(E != NULL);
    size_t * P = malloc(N*sizeof(double)); // Keep also bead numbers
    assert(P!=NULL);

    /* Move data into new structure */
    for(size_t kk = 0; kk<N; kk++)
    {
        size_t h = hash(nDiv, D+3*kk);

        size_t writepos = 3*C[h];
        C[h]++;
        P[writepos/3] = kk;
        for(int idx = 0; idx<3; idx++)
        {
            assert(writepos<N*3);
            E[writepos++] = D[3*kk+idx];
        }
    }

    /*
     * now we can hash(X) a point X,
     * using B to see where the bucket is
     * in E
     * */

    // Get some job done!
    double err = 0;
    for(size_t kk = 0; kk<N; kk++)
        /* Looping over E in order to avoid self-matches
         * and duplicates */
    {
        double deps = d;
        size_t ha_min = hash_coord(nDiv, E[3*kk]-deps);
        size_t ha_max = hash_coord(nDiv, E[3*kk]+deps);
        size_t hb_min = hash_coord(nDiv, E[3*kk+1]-deps);
        size_t hb_max = hash_coord(nDiv, E[3*kk+1]+deps);
        size_t hc_min = hash_coord(nDiv, E[3*kk+2]-deps);
        size_t hc_max = hash_coord(nDiv, E[3*kk+2]+deps);

        for(size_t cc = hc_min; cc <= hc_max; cc++)
            for(size_t bb = hb_min; bb <= hb_max; bb++)
                for(size_t aa = ha_min; aa <= ha_max; aa++)
                {

                    size_t hash =
                        aa +
                        bb*nDiv +
                        cc*pow(nDiv,2);
                    for(size_t pp = B[hash]; pp<B[hash+1]; pp++)
                    {
                        if(pp>kk) {
                            double dist2 = eudist3sq(E+3*pp, E+3*kk);
                            if(dist2<d2)
                            {

                                // Retrieve original positions
                                size_t KK = P[kk];
                                size_t LL = P[pp];

                                double di = sqrt(dist2);

                                double did = (di - d)/di;

                                if( !isfinite(did))
                                {
                                    did = 0;
                                }

                                for(int idx = 0; idx<3; idx++)
                                {
                                    G[3*KK+idx] += kVol*2*(E[3*kk+idx] - E[3*pp+idx])*did;
                                    G[3*LL+idx] -= kVol*2*(E[3*kk+idx] - E[3*pp+idx])*did;
                                }

                            }
                        }
                    }
                }
    }


    free(B);
    free(C);
    free(E);
    free(S);
    free(P);

    return err;
}


double err3(const double * restrict X,
            const size_t nX,
            const double * restrict R,
            const uint32_t * restrict P,
            const mflock_func_t * restrict C )
{
    /* Compared to err2 this is an alternative version with a list of pairs in contact instead of A */

    double XT[3]; // Temporary storage
    XT[0] = 0; XT[1] = 0; XT[2] = 0;

    /* -- Wanted radii -- */
    double errRad = 0;
    if(C->E == NULL) // spherical domain
    {
        if( (R != NULL) & (C->kRad > 0))
        {
            for(size_t kk = 0; kk<nX; kk++)
            {
                if(isfinite(R[kk]) == 1) // TODO: create list of points to use to avoid if-branching
                {
                    const double r = norm3(X+kk*3);
                    errRad += pow(r-R[kk], 2);
                }
            }
        }
    }
    if(C->E != NULL) // elliptical domain
    {
        if( (R != NULL) & (C->kRad > 0))
        {
            for(size_t kk = 0; kk<nX; kk++)
            {
                if(isfinite(R[kk]) == 1)
                {
                    const double r = sqrt(elli_getScale2(C->E, X+3*kk));
                    errRad += pow(r-R[kk], 2);
                }
            }
        }
    }

    /* -- Keep inside domain -- */
    double errSph = 0;

    if(C->E == NULL)
    { // Inside sphere
        double ds = (1-C->r0);
        double ds2 = pow(ds, 2);
        for(size_t kk = 0 ; kk<nX; kk++)
        {
            const double r2 = norm32(X+kk*3); // norm^2
            if( r2 > ds2)
            {
                double r = sqrt(r2);
                errSph += pow(r - ds, 2);
            }
        }
    }
    if(C->E != NULL) // Inside ellipsoid
    {
        for(size_t kk = 0 ; kk<nX; kk++)
        {
            const double d = elli_getScale2(C->Es, X+3*kk);
            if( d >= 1)
            {
                double r = elli_gdistL(C->E, X+3*kk, XT);
                if(r + C->r0 > 0)
                {
                    errSph += pow(r+C->r0, 2);
                }
            }
        }
    }


    /* Wanted contacts/interactions */
    double errInt = 0;
    for(size_t pp = 0; pp < C->nIPairs; pp++)
    {
        size_t kk = P[pp*2];
        size_t ll = P[pp*2+1];
        {
            double d = eudist3(X+3*kk, X+3*ll);
            // TODO: min/max to avoid branch?
            if(d > C->dInteraction)
            {
                errInt += pow(d - C->dInteraction, 2);
            }
        }
    }

    // Repulsion
    double errVol = 0;

    errVol = errRepulsion(X, nX, 2*C->r0);

    return C->kInt*errInt + C->kVol*errVol + C->kDom*errSph + C->kRad*errRad;
}

double err2(double * X, size_t nX, double * R, uint32_t * P, mflock_func_t * C )
{
    /* Alternative version with a list of pairs in contact instead of A */

    // Wanted radii
    double errRad = 0;
    if( (R != NULL) & (C->kRad > 0))
    {
        for(size_t kk = 0; kk<nX; kk++)
        {
            if(isfinite(R[kk]) == 1)
            {
                const double r = norm3(X+kk*3);
                errRad += pow(r-R[kk], 2);
            }
        }
        //  printf("errRad: %f\n", errRad);
    }


    // Keep inside sphere
    double errSph = 0;
    for(size_t kk = 0 ; kk<nX; kk++)
    {
        const double r = norm3(X+kk*3);
        if( r > 1-C->r0)
        {
            errSph += pow(r + C->r0 -1, 2);
        }
    }

    // Wanted contacts/interactions
    double errInt = 0;
    for(size_t pp = 0; pp < C->nIPairs; pp++)
    {
        size_t kk = P[pp*2];
        size_t ll = P[pp*2+1];
        {
            double d = eudist3(X+3*kk, X+3*ll);
            if(d > C->dInteraction)
            {
                errInt += pow(d - C->dInteraction, 2);
            }
        }
    }

    // Repulsion
    double errVol = 0;
    for(size_t kk = 0; kk < nX; kk++)
    {
        for(size_t ll = kk+1; ll < nX; ll++)
        {
            if(fabs(X[3*kk] - X[3*ll]) < 2*C->r0) // This line doubles the speed!
            {
                double d = eudist3(X+3*kk, X+3*ll);
                if( d < 2*C->r0)
                {
                    errVol += pow(d - 2*C->r0, 2);
                }
            }
        }
    }

    return C->kInt*errInt + C->kVol*errVol + C->kDom*errSph + C->kRad*errRad;
}


double err(double * X, size_t nX, double * R, uint8_t * A, mflock_func_t * C )
{

    // Wanted radii
    double errRad = 0;
    if( (R != NULL) & (C->kRad > 0))
    {
        for(size_t kk = 0; kk<nX; kk++)
        {
            if(isfinite(R[kk]) == 1)
            {
                const double r = norm3(X+kk*3);
                errRad += pow(r-R[kk], 2);
            }
        }
        //    printf("errRad: %f\n", errRad);
    }

    // Keep inside sphere
    double errSph = 0;
    for(size_t kk = 0 ; kk<nX; kk++)
    {
        const double r = norm3(X+kk*3);
        if( r > 1-C->r0)
        {
            errSph += pow(r + C->r0 -1, 2);
        }
    }

    // Wanted contacts/interactions
    double errInt = 0;
    for(size_t kk = 0; kk < nX; kk++)
    {
        for(size_t ll = kk+1; ll < nX; ll++)
        {
            if(A[kk + nX*ll] == 1)
            {
                double d = eudist3(X+3*kk, X+3*ll);
                if(d > C->dInteraction)
                {
                    errInt += pow(d - C->dInteraction, 2);
                }
            }
        }
    }

    // Repulsion
    double errVol = 0;
    for(size_t kk = 0; kk < nX; kk++)
    {
        for(size_t ll = kk+1; ll < nX; ll++)
        {
            double d = eudist3(X+3*kk, X+3*ll);
            if( d < 2*C->r0)
            {
                errVol += pow(d - 2*C->r0, 2);
            }
        }
    }

    return C->kInt*errInt + C->kVol*errVol + C->kDom*errSph + C->kRad*errRad;
}

void grad(double * X, size_t nX, double * R, uint8_t * A, double * G, mflock_func_t * C)
{
    // Reset G
    for(size_t kk = 0; kk<nX*3; kk++)
        G[kk] = 0;

    // Radial positioning
    if(C->kRad > 0)
    {
        for(size_t kk = 0; kk<nX; kk++)
        {
            if(isfinite(R[kk]) == 1)
            {
                double r = norm3(X+kk*3);
                double re = 0;
                if(r > 0)
                    re = 2*1/r*(r-R[kk]);
                for(int idx = 0; idx<3; idx++)
                {
                    G[3*kk+idx] += C->kRad*X[kk*3+idx]*re;
                }
            }
        }
    }

    // Keep inside sphere
    for(size_t kk = 0; kk<nX; kk++)
    {
        double r = norm3(X+kk*3);
        if(r > 1-C->r0)
        {
            double re = 2 / r * (r-(1-C->r0));
            for(int idx = 0; idx<3; idx++)
            {
                G[3*kk+idx] += C->kDom*X[kk*3+idx]*re;
            }
        }
    }

    // Wanted interactions
    for(size_t kk = 0; kk < nX; kk++)
    {
        for(size_t ll = kk+1; ll < nX; ll++)
        {
            if(A[kk + nX*ll] == 1)
            {
                double d = eudist3(X+3*kk, X+3*ll);
                if(d > C->dInteraction)
                {
                    for(int idx = 0; idx<3; idx++)
                    {
                        G[3*kk+idx] += C->kInt*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - C->dInteraction);
                        G[3*ll+idx] -= C->kInt*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - C->dInteraction);
                    }
                }
            }
        }
    }

    // Repulsion
    for(size_t kk = 0; kk < nX; kk++)
    {
        for(size_t ll = kk+1; ll < nX; ll++)
        {
            double d = eudist3(X+3*kk, X+3*ll);
            if( d < 2*C->r0)
            {
                for(int idx = 0; idx<3; idx++)
                {
                    G[3*kk+idx] += C->kVol*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - 2*C->r0);
                    G[3*ll+idx] -= C->kVol*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - 2*C->r0);
                }
            }
        }
    }

    return;
}

void grad2(double * X, size_t nX, double * R, uint32_t * I, double * G, mflock_func_t * C)
{
    // Reset G
    for(size_t kk = 0; kk<nX*3; kk++)
        G[kk] = 0;

    // Radial positioning
    if(C->kRad > 0)
    {
        for(size_t kk = 0; kk<nX; kk++)
        {
            if(isfinite(G[kk]) == 1)
            {
                double r = norm3(X+kk*3);
                double re = 0;
                if(r > 0)
                    re = 2*1/r*(r-R[kk]);
                for(int idx = 0; idx<3; idx++)
                {
                    G[3*kk+idx] += C->kRad*X[kk*3+idx]*re;
                }
            }
        }
    }



    // Keep inside sphere
    for(size_t kk = 0; kk<nX; kk++)
    {
        double r = norm3(X+kk*3);
        if(r > 1-C->r0)
        {
            double re = 2 / r * (r-(1-C->r0));
            for(int idx = 0; idx<3; idx++)
            {
                G[3*kk+idx] += C->kDom*X[kk*3+idx]*re;
            }
        }
    }

    // Wanted interactions
    for(size_t pp = 0; pp < C->nIPairs; pp++)
    {
        size_t kk = I[pp*2];
        size_t ll = I[pp*2+1];

        double d = eudist3(X+3*kk, X+3*ll);
        if(d > C->dInteraction)
        {
            for(int idx = 0; idx<3; idx++)
            {
                G[3*kk+idx] += C->kInt*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - C->dInteraction);
                G[3*ll+idx] -= C->kInt*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - C->dInteraction);
            }
        }
    }

    // Repulsion
    for(size_t kk = 0; kk < nX; kk++)
    {
        for(size_t ll = kk+1; ll < nX; ll++)
        {
            double d = eudist3(X+3*kk, X+3*ll);
            if( d < 2*C->r0)
            {
                for(int idx = 0; idx<3; idx++)
                {
                    G[3*kk+idx] += C->kVol*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - 2*C->r0);
                    G[3*ll+idx] -= C->kVol*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - 2*C->r0);
                }
            }
        }
    }

    return;
}

void grad3(const double * restrict X,
           const size_t nX,
           const double * restrict R,
           const uint32_t * restrict I,
           double * restrict G,
           const mflock_func_t * restrict C)
{

    double XT[3]; XT[0] = 0; XT[1] = 0; XT[2] = 0;

    memset(G, 0, nX*3*sizeof(double));

    /* Radial positioning */
    if(C->E == NULL) // Spherical domain
    {
        if(C->kRad > 0)
        {
            for(size_t kk = 0; kk<nX; kk++)
            {
                if(isfinite(R[kk]) == 1)
                {
                    double r = norm3(X+kk*3);
                    double re = 0;
                    // if(r > 0)
                    re = 2*1/r*(r-R[kk]);
                    for(int idx = 0; idx<3; idx++)
                    {
                        G[3*kk+idx] += C->kRad*X[kk*3+idx]*re;
                    }
                }
            }
        }
    }

    if(C->E != NULL) // Ellipsoidal domain
    {
        if(C->kRad > 0)
        {
            double EF[3];
            EF[0] = 1.0/pow(C->E->a, 2);
            EF[1] = 1.0/pow(C->E->b, 2);
            EF[2] = 1.0/pow(C->E->c, 2);

            for(size_t kk = 0; kk<nX; kk++)
            {
                if(isfinite(R[kk]) == 1)
                {
                    double r = sqrt(elli_getScale2(C->E, X+kk*3));
                    double re = 0;
                    // if(r > 0)
                    re = 2*(r-R[kk])/r;
                    for(int idx = 0; idx<3; idx++)
                    {
                        G[3*kk+idx] += C->kRad*X[kk*3+idx]*re*EF[idx];
                    }
                }
            }
        }
    }

    /* --- Keep inside domain --- */

    if(C->E == NULL) // Spherical domain
    {
        for(size_t kk = 0; kk<nX; kk++)
        {
            double r = norm3(X+kk*3);
            if(r > 1-C->r0)
            {
                double re = 2 / r * (r-(1-C->r0));
                for(int idx = 0; idx<3; idx++)
                {
                    G[3*kk+idx] += C->kDom*X[kk*3+idx]*re;
                }
            }
        }
    }

    if(C->E != NULL) // Ellipsoidal domain
    {
        for(size_t kk = 0; kk<nX; kk++)
        {
            const double d = elli_getScale2(C->Es, X+3*kk);
            if( d >= 1)
            {
                double n[3];
                double r = elli_gdistL(C->E, X+3*kk, XT);

                //elli_normal(C->E, XT, n);
                vec3_normalize(n);

                if(r + C->r0 > 0)
                {
                    double re = 2 / r * (r+C->r0);
                    for(int idx = 0; idx<3; idx++)
                    {
                        G[3*kk+idx] += C->kDom*(X[kk*3+idx]-XT[idx])*re;

                    }
                }
            }
        }
    }

    // Wanted interactions
    for(size_t pp = 0; pp < C->nIPairs; pp++)
    {
        size_t kk = I[pp*2];
        size_t ll = I[pp*2+1];

        double d = eudist3(X+3*kk, X+3*ll);
        assert(kk != ll);
#ifndef NDEBUG
        if( !(d>0) )
        {
            printf("Strange distance between interacting points!\n");
            printf("@%p : %f %f %f\n", (void*) (X+3*kk), X[3*kk], X[3*kk+1], X[3*kk+2]);
            printf("@%p : %f %f %f\n", (void*) (X+3*ll), X[3*ll], X[3*ll+1], X[3*ll+2]);
            exit(1);
        }
#endif
        if(d > C->dInteraction && d > 1e-6)
        {
            for(int idx = 0; idx<3; idx++)
            {
                G[3*kk+idx] += C->kInt*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - C->dInteraction);
                G[3*ll+idx] -= C->kInt*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - C->dInteraction);
            }
        }
    }

    // Repulsion
    gradRepulsion(X, G, nX, 2*C->r0, C->kVol);

    return;
}


void grad4(double * restrict X,
           const size_t nX,
           double * restrict R,
           uint32_t * restrict I,
           double * restrict G,
           const mflock_func_t * restrict C)
{
    // Reset G
    for(size_t kk = 0; kk<nX*3; kk++)
        G[kk] = 0;

    // Radial positioning
    if(C->kRad > 0)
    {
        for(size_t kk = 0; kk<nX; kk++)
        {
            if(isfinite(R[kk]) == 1)
            {


                double r = norm3(X+kk*3);
                double re = 0;
                // if(r > 0)
                re = 2*1/r*(r-R[kk]);
                for(int idx = 0; idx<3; idx++)
                {
                    G[3*kk+idx] += C->kRad*X[kk*3+idx]*re;
                }
            }
        }
    }

    // Keep inside sphere
    const double rmax = 1-C->r0;
    for(size_t kk = 0; kk<nX; kk++)
    {
        double r = norm3(X+kk*3);
        if(r > rmax)
        {
            double re = 2 / r * (r-(1-C->r0));
            for(int idx = 0; idx<3; idx++)
            {
                G[3*kk+idx] += C->kDom*X[kk*3+idx]*re;
            }
        }
    }

    // Wanted interactions
    const double dInteraction = C->dInteraction;
    const double dInteraction2 = pow(dInteraction, 2);

    for(size_t pp = 0; pp < C->nIPairs; pp++)
    {
        size_t kk = I[pp*2];
        size_t ll = I[pp*2+1];

        double d2 = eudist3sq(X+3*kk, X+3*ll);
        assert(kk != ll);
        assert(d2>0);
        if(d2 > dInteraction2)
        {
            double d= sqrt(d2);

            for(int idx = 0; idx<3; idx++)
            {
                G[3*kk+idx] += C->kInt*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - dInteraction);
                G[3*ll+idx] -= C->kInt*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - dInteraction);
            }
        }
    }

    // Repulsion
    gradRepulsion(X, G, nX, 2*C->r0, C->kVol);

    return;
}


void bead_wells_gradient(const mflock_func_t * restrict fconf,
                         const double * restrict W,
                         const size_t nW,
                         const double * restrict X,
                         double * restrict G)
{
    double sigma = fconf->r0;
    double K1 = 1.0 / sigma; // 1.0 / pow(sigma, 3) / sqrt(2.0*M_PI);
    K1 *= fconf->kBeadWell;
    double K2 = -0.5/pow(sigma, 2);

    for(size_t kk = 0; kk < nW; kk++)
    {
        size_t bead = (size_t )W[4*kk];
        const double * well = W + 4*kk + 1;
        const double * pos = X + 3*bead;

        double r2 = pow(pos[0]-well[0], 2) +
            pow(pos[1]-well[1], 2) +
            pow(pos[2]-well[2], 2);

        double K3 = K1*exp(K2*r2);

        for(size_t ii = 0; ii< 3; ii++)
        {
            double delta_i = pos[ii]-well[ii];
            G[3*bead+ii] += delta_i*K3;
        }
    }
    return;
}

double
bead_wells_error(const mflock_func_t * restrict fconf,
                 const double * restrict W,
                 const size_t nW,
                 const double * restrict X)
{

    double E = 0;
    double sigma = fconf->r0;

    const double c0 = fconf->kBeadWell; // / sigma / sqrt(2.0*M_PI);

    double K1 = 1.0; // 1.0 / sigma / sqrt(2.0*M_PI);

    K1 *= fconf->kBeadWell;
    double K2 = -0.5/pow(sigma, 2);

    for(size_t kk = 0; kk < nW; kk++)
    {
        size_t bead = (size_t )W[4*kk];
        const double * well = W + 4*kk + 1;
        const double * pos = X + 3*bead;

        double r2 = pow(pos[0]-well[0], 2) +
            pow(pos[1]-well[1], 2) +
            pow(pos[2]-well[2], 2);

        E += (c0 - K1*exp(K2*r2));
    }
    return E;
}
