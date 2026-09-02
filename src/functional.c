#include "functional.h"


static double f_interaction(double r)
{
    return pow(r, 2.0);
}

/* Squared L2 norm of a 3-vector */
static double norm32(const double * restrict X)
{
    double n = 0;
    for(size_t kk = 0; kk<3; kk++)
    {
        n+=pow(X[kk], 2);
    }
    return n;
}

static double norm3(const double * restrict X)
{
    return sqrt(norm32(X));
}

static void vec3_normalize(double * restrict X)
{
    double n = norm3(X);
    X[0] /= n; X[1] /= n; X[2] /= n;
}


static double eudist3sq(const double * A, const double * B)
{
    /* SQUARED Euclidean distance between two 3D-vectors */
    return pow(A[0]-B[0], 2) + pow(A[1]-B[1], 2) + pow(A[2]-B[2],2);
}

static double eudist3(const double * A, const double * B)
{
    /* Euclidean distance between two 3D-vectors */
    return sqrt(eudist3sq(A, B));
}

typedef struct {
    const double * X;
    double * G;
    double radius;
    double kVol;
} repulsion_gradient_params;

typedef struct {
    const double * X;
    double error;
    double radius;
    double kVol;
} repulsion_error_params;

void repulsion_gradient_cb(u32 u, u32 v, double d2, void * data)
{
    repulsion_gradient_params * params = (repulsion_gradient_params *) data;
    if( d2 < 1e-8 ) { // avoid division by zero
        return;
    }
    double d = sqrt(d2);

    double * G = params->G;
    const double * X = params->X;
    double kVol = params->kVol;

    double did = (d - params->radius)/d;
    double f = 2*kVol*did;

    for(i64 ii = 0; ii < 3; ii++){
        G[3*u + ii] += f*(X[3*u + ii] - X[3*v + ii]);
        G[3*v + ii] -= f*(X[3*u + ii] - X[3*v + ii]);
    }
}

void repulsion_error_cb(__attribute__((unused)) u32 u,
                        __attribute__((unused)) u32 v,
                        double d2, void * data)
{
    repulsion_error_params * params = (repulsion_error_params *) data;
    if( d2 < 1e-8 ) { // avoid division by zero
        return;
    }
    double d = sqrt(d2);
    params->error += pow(d - params->radius, 2);
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
        if( (R != NULL) && (C->kRad > 0))
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
        if( (R != NULL) && (C->kRad > 0))
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
            if(d > C->dInteraction)
            {
                //errInt += pow(d - C->dInteraction, 2);
                errInt += f_interaction(d - C->dInteraction);
            }
        }
    }

    // Repulsion
    collide3_info cinfo = {};
    repulsion_error_params cparams = {};
    cparams.error = 0;
    cparams.radius = 2.0*C->r0;
    cparams.kVol = C->kVol;
    cparams.X = X;
    collide3_f64(X, nX,
                 2*C->r0, &cinfo,
                 repulsion_error_cb,
                 &cparams);
    double errVol = cparams.error;

    return C->kInt*errInt + C->kVol*errVol + C->kDom*errSph + C->kRad*errRad;
}

double err2(double * X, size_t nX, double * R, uint32_t * P, mflock_func_t * C )
{
    /* Alternative version with a list of pairs in contact instead of A */

    // Wanted radii
    double errRad = 0;
    if( (R != NULL) && (C->kRad > 0))
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
                //errInt += pow(d - C->dInteraction, 2);
                errInt += f_interaction(d - C->dInteraction);
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
    if( (R != NULL) && (C->kRad > 0))
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
                    //errInt += pow(d - C->dInteraction, 2);
                    errInt += f_interaction(d - C->dInteraction);
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


void
grad3(const double * restrict X,
      const size_t nX,
      const double * restrict R,
      const uint32_t * restrict I,
      uint8_t * restrict active_pair,
      const uint32_t * restrict backbone,
      const size_t n_backbone,
      double * restrict G,
      const mflock_func_t * restrict C)
{

    double XT[3]; XT[0] = 0; XT[1] = 0; XT[2] = 0;

    memset(G, 0, nX*3*sizeof(double));

    /* Radial positioning */
    if( (C->kRad > 0) && (R != NULL)) {
        if(C->geometry == MFLOCK_ELLIPSOID)
        {
            printf("Warning: Using radial constrains with ellipsoidal geometry"
                   "is not implemented yet. Will give you weird results\n");
            // Obviously we need to convert the radial values to something else
            // in this case.
        }
        for(size_t kk = 0; kk<nX; kk++) {
            if(isfinite(R[kk]) == 1) {
                double r = norm3(X+kk*3);
                double re = 0;
                // if(r > 0)
                re = 2*1/r*(r-R[kk]);
                for(int idx = 0; idx<3; idx++) {
                    G[3*kk+idx] += C->kRad*X[kk*3+idx]*re;
                }
            }
        }
    }


    if(C->geometry == MFLOCK_ELLIPSOID)
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
    if(C->geometry == MFLOCK_SPHERE)
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

    // Top and bottom restrictions
    if(C->top_plane < 1.0) {

        for(size_t kk = 0; kk<nX; kk++) {

            if( X[3*kk + 2] > C->top_plane ) {
                G[3*kk + 2] -= C->kDom*(C->top_plane - X[kk*3 + 2]);
            }
        }
    }

    if(C->bottom_plane > -1.0) {

        for(size_t kk = 0; kk<nX; kk++) {

            if( X[3*kk + 2] < C->bottom_plane ) {
                G[3*kk + 2] -= C->kDom*(C->bottom_plane - X[kk*3 + 2]);
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
                double n[3] = {0};
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
    if(active_pair == NULL)
    {
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
    } else {
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

            if(active_pair[pp] == 0){
                if(d < 3.0*C->r0) {
                    active_pair[pp] = 1;
                }
            }

            if(active_pair[pp] == 1) {
            if(d > C->dInteraction && d > 1e-6)
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
    // backbone
    for(size_t pp = 0; pp < n_backbone; pp++)
    {
        size_t kk = backbone[pp*2];
        size_t ll = backbone[pp*2+1];

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
                G[3*kk+idx] += C->kBackbone*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - C->dInteraction);
                G[3*ll+idx] -= C->kBackbone*2*(X[3*kk+idx] - X[3*ll+idx])/d*(d - C->dInteraction);
            }
        }
    }

    // Repulsion
    collide3_info cinfo = {};
    repulsion_gradient_params cparams;
    cparams.G = G;
    cparams.radius = 2.0*C->r0;
    cparams.kVol = C->kVol;
    cparams.X = X;

    collide3_f64(X, nX,
                 2*C->r0, &cinfo,
                 repulsion_gradient_cb,
                 &cparams);

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
    collide3_info cinfo = {};
    repulsion_gradient_params cparams;
    cparams.G = G;
    cparams.radius = 2.0*C->r0;
    cparams.kVol = C->kVol;
    cparams.X = X;

    collide3_f64(X, nX,
                 2*C->r0, &cinfo,
                 repulsion_gradient_cb,
                 &cparams);

    return;
}


void
bead_wells_gradient(const mflock_func_t * restrict fconf,
                    const size_t n_bead,
                    const wpos * restrict W,
                    const size_t nW,
                    const double * restrict X,
                    double * restrict G)
{
    double sigma = fconf->r0;
    double K1 = 1.0 / sigma; // 1.0 / pow(sigma, 3) / sqrt(2.0*M_PI);
    K1 *= fconf->kBeadWell;
    double K2 = -0.5/pow(sigma, 2);

    if(fconf->diploid == 0)
    {
        for(size_t kk = 0; kk < nW; kk++)
        {
            const wpos well = W[kk];
            const double * pos = X + 3*well.bead_idx;

            double r2 = pow(pos[0]-well.X[0], 2) +
                pow(pos[1]-well.X[1], 2) +
                pow(pos[2]-well.X[2], 2);

            double K3 = K1*exp(K2*r2);
            if(0){
                printf("%zu: (%f, %f, %f) -> (%f, %f, %f) (%f)\n",
                       well.bead_idx,
                       pos[0], pos[1], pos[2],
                       well.P.x, well.P.y, well.P.z,
                       K3);
            }
            for(size_t ii = 0; ii< 3; ii++)
            {
                double delta_i = pos[ii]-well.X[ii];
                G[3*well.bead_idx+ii] += delta_i*K3;
            }
        }
    } else {
        assert(fconf->diploid == 1);
        for(size_t kk = 0; kk < nW; kk++)
        {
            const wpos well = W[kk];

            const double * posA = X + 3*well.bead_idx;

            double r2A = pow(posA[0]-well.X[0], 2) +
                pow(posA[1]-well.X[1], 2) +
                pow(posA[2]-well.X[2], 2);

            const double * posB = X + 3*(well.bead_idx + n_bead/2);

            double r2B = pow(posB[0]-well.X[0], 2) +
                pow(posB[1]-well.X[1], 2) +
                pow(posB[2]-well.X[2], 2);

            const double * pos = posA;
            double r2 = r2A;

            // index of bead that should receive the gradient
            size_t idx = well.bead_idx;
            if(r2B < r2A)
            {
                r2 = r2B;
                pos = posB;
                idx += n_bead/2;
            }

            double K3 = K1*exp(K2*r2);
            for(size_t ii = 0; ii< 3; ii++)
            {
                double delta_i = pos[ii]-well.X[ii];
                G[3*idx+ii] += delta_i*K3;
            }
        }
    }
    return;
}

double
bead_wells_error(const mflock_func_t * restrict fconf,
                 const size_t n_bead,
                 const wpos * restrict W,
                 const size_t nW,
                 const double * restrict X)
{

    double E = 0;
    double sigma = fconf->r0;

    const double c0 = fconf->kBeadWell; // / sigma / sqrt(2.0*M_PI);

    double K1 = 1.0; // 1.0 / sigma / sqrt(2.0*M_PI);
    K1 *= fconf->kBeadWell;
    double K2 = -0.5/pow(sigma, 2);

    if(fconf->diploid == 0)
    {

        for(size_t kk = 0; kk < nW; kk++)
        {
            wpos well = W[kk];
            const double * pos = X + 3*well.bead_idx;

            double r2 = pow(pos[0]-well.X[0], 2) +
                pow(pos[1]-well.X[0], 2) +
                pow(pos[2]-well.X[0], 2);

            E += (c0 - K1*exp(K2*r2));
        }
    } else {
        assert(fconf->diploid == 1);
        for(size_t kk = 0; kk < nW; kk++)
        {
            wpos well = W[kk];
            const double * posA = X + 3*well.bead_idx;
            double r2A = pow(posA[0]-well.X[0], 2) +
                pow(posA[1]-well.X[0], 2) +
                pow(posA[2]-well.X[0], 2);
            const double * posB = X + 3*(well.bead_idx+n_bead/2);
            double r2B = pow(posB[0]-well.X[0], 2) +
                pow(posB[1]-well.X[0], 2) +
                pow(posB[2]-well.X[0], 2);
            double r2 = r2A;
            r2 > r2B ? r2 = r2A : 0;
            E += (c0 - K1*exp(K2*r2));
        }
    }
    return E;
}
