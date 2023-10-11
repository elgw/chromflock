#include "mflock.h"


static volatile int run = 1;

#define INLINED static inline __attribute__((always_inline))

static void luaerror (lua_State *L, const char *fmt, ...) {
    /* show lua error string */
    va_list argp;
    va_start(argp, fmt);
    vfprintf(stderr, fmt, argp);
    printf("\n");
    va_end(argp);
    lua_close(L);
    exit(EXIT_FAILURE);
}

static double getglobfloat (lua_State *L, const char *var) {
    /* get global variable by name */
    int isnum;
    double result;
    lua_getglobal(L, var);
    result = (double)lua_tonumberx(L, -1, &isnum);
    if (!isnum)
        luaerror(L, "'%s' should be a number\n", var);
    lua_pop(L, 1);
    /* remove result from the stack */
    return result;
}

static int getglobint (lua_State *L, const char *var) {
    /* get global variable by name */
    int isnum, result;
    lua_getglobal(L, var);
    result = (int)lua_tointegerx(L, -1, &isnum);
    if (!isnum)
        luaerror(L, "'%s' should be a number\n", var);
    lua_pop(L, 1);
    /* remove result from the stack */
    return result;
}

static void stoprun(int ignore)
{
    run = 0 + 0*ignore;
}

static double clockdiff(struct timespec* start, struct timespec * finish)
{
    double elapsed = (finish->tv_sec - start->tv_sec);
    elapsed += (finish->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}

static double norm3(const double * restrict X)
{
    double n = 0;
    for(size_t kk = 0; kk<3; kk++)
        n+=pow(X[kk], 2);
    return sqrt(n);
}

INLINED double dmax(double a, double b)
{
    if(a>b)
        return a;
    return b;
}

INLINED double dmin(double a, double b)
{
    if(a<b)
        return a;
    return b;
}

typedef struct {
    double * X;
    optparam * p;
} solve_t_args;

void radpos(double * X, optparam * p, FILE * f)
{
    /* Show radial positions of each chromosome, using centre of mass
     * an alternative would be to show mean radius. */

    if(p->L == NULL)
    {
        printf("No L!\n");
        return;
    }
    // mX mean X positions
    size_t size_mX = 3*256*sizeof(double);
    double * mX = malloc(size_mX);
    assert(mX != NULL);
    memset(mX, 0, size_mX);
    // Number of points per label
    size_t * nL = malloc(256*sizeof(size_t));
    assert(nL != NULL);
    memset(nL, 0, 256*sizeof(size_t));

    // Get centroid for each chromosome/label
    for(size_t kk = 0 ; kk < p->N ; kk++)
    {
        size_t label = p->L[kk];
        nL[label]++;
        for(size_t idx = 0; idx < 3 ; idx ++)
        {
            mX[3*label+idx] += X[3*kk+idx];
        }
    }

    // Print the results
    fprintf(f, "chr, com_x, com_y, com_z, r\n");
    for(size_t ll = 0; ll < 256; ll++)
    {
        if(nL[ll] > 0)
        {
            for(size_t idx = 0; idx<3; idx++)
            {
                mX[3*ll+idx]/=nL[ll];
            }
            fprintf(f, "%3zu,% 2.3f, % 2.3f, % 2.3f, %.3f\n", ll,
                    mX[3*ll],
                    mX[3*ll+1],
                    mX[3*ll+2],
                    norm3(mX+3*ll));
        }
    }

    free(mX);
    free(nL);
    return;
}

/* when usleep is called from the lua script
 * This is useful when visualizations are enabled
 */
static int usleep_lua(lua_State * L)
{
    /* get number of arguments */
    if(lua_gettop(L) != 1)
    {
        lua_pushstring(L, "Incorrect number of arguments to 'usleep'");
        lua_error(L);
        return 0;
    }

    if (!lua_isnumber(L, 1))
    {
        lua_pushstring(L, "Incorrect argument to 'usleep'");
        lua_error(L);
        return 0;
    }

    double utime = lua_tonumber(L, 1);
    usleep((size_t) utime);

    /* return the number of results */
    return 0;

}

INLINED double eudist3p2(const double * A, const double * B)
{
    /* Euclidean distance between two 3D-vectors */
    return pow(A[0]-B[0], 2) + pow(A[1]-B[1], 2) + pow(A[2]-B[2],2);
}

INLINED double eudist3(const double * A, const double * B)
{
    /* Euclidean distance between two 3D-vectors */
    return sqrt( pow(A[0]-B[0], 2) + pow(A[1]-B[1], 2) + pow(A[2]-B[2], 2));
}


/** @brief Per Chr Centre of mass compression
 *
 *Apply a force that attracts each bead to the centre of mass of it's
 * chromosome This slows down the computations quite much because the
 * beads get close to each other so the collision detection gets more to do.
 *
 * TODO: Unnecessary to allocate/free things here and to count the
 * number of beads per chromosome.
 */
void comforce(optparam * restrict p,
              const double * restrict X,
              double * restrict G)
{
    if(p->L == NULL)
    {
        printf("No L, can't compute com-force!\n");
        return;
    }

    /* mX mean X positions */
    size_t size_mX = 3*256*sizeof(double);
    double * mX = malloc(size_mX);
    memset(mX, 0, size_mX);

    /* Number of points per label */
    size_t * nL = malloc(256*sizeof(size_t));
    memset(nL, 0, 256*sizeof(size_t));

    /* Get centroid for each chromosome/label */
    for(size_t kk = 0 ; kk < p->N ; kk++)
    {
        size_t label = p->L[kk];
        nL[label]++;
        for(size_t idx = 0; idx < 3 ; idx ++)
        {
            mX[3*label+idx] += X[3*kk+idx];
        }
    }

    /* Normalize by the number of points */
    for(size_t ll = 0 ; ll < 256 ; ll++)
    {
        double nPoints = nL[ll];
        if(nPoints > 0)
        {
            for(size_t idx = 0; idx < 3 ; idx ++)
            {
                mX[3*ll+idx] /= nPoints;
            }
        }
    }

    /* Attract to centroid */
    for(size_t kk = 0; kk < p->N; kk++)
    {
        uint8_t label = p->L[kk];
        double dist2 = eudist3p2(X+3*kk, mX+3*label);
        if(dist2 > 0)
        {
            for(size_t idx = 0; idx<3; idx++)
            {
                G[3*kk+idx] += p->compress*(X[3*kk+idx]-mX[3*label + idx])*dist2;
            }
        }
    }

    free(nL);
    free(mX);
    return;
}

/**
 * @breif The beads dynamics main loop
 * @param X the coordinates of each bead
 * @param p the settings
 * @param Fb: the Brownian force
 */
int dynamic(double * restrict X,
            optparam * restrict p,
            double Fb)
{

    normal_setup();  /* set up fast prng */

    /* Prepare settings */
    size_t maxiter = p->maxiter;
    conf fconf;
    fconf.r0 = p->r0;
    fconf.kVol = p->kVol;
    fconf.kDom = p->kDom;
    fconf.kInt = p->kInt; /* Function of current iteration, see below */
    fconf.kRad = p->kRad;
    fconf.dInteraction = 2.1*fconf.r0; /* Function of current iteration, see below */
    fconf.E = NULL;
    fconf.Es = NULL;
    if(p->E != NULL)
    {
        /* add ellipse parameters otherwise sphere domain */
        fconf.E = p->E;
        fconf.Es = elli_new(
                            p->E->a - fconf.r0,
                            p->E->b - fconf.r0,
                            p->E->c - fconf.r0);
    }

    fconf.nIPairs = p->NI; /* Only for err2 */

    /* Initialization */
    struct timespec ta, tb;
    clock_gettime(CLOCK_MONOTONIC, &ta);

    /* Gradient */
    double * restrict g = malloc(p->N*3*sizeof(double));
    memset(g, 0, p->N*3*sizeof(double)); /* zero initial velocity */
    // Velocity
    double * restrict v = malloc(p->N*3*sizeof(double));
    memset(v, 0, p->N*3*sizeof(double)); /* zero initial velocity */
    double * restrict Xm = malloc(3*p->N*sizeof(double));
    /* Set last position to be the current position, i.e. initial velocity will be 0 */
    for(size_t pp = 0; pp< 3*p->N ; pp++)
    { Xm[pp] = X[pp]; }

    double error = 9e99;
    double gnorm = 9e99;
    double dt = 0.15;
    double damp = 0.5; /* dampening (.8) faster convergence than .4 */
    int showStep = 0;
    double Frnd = .10; /* Don't change this, change Fb instead */

    lua_State *L = NULL;

    if(p->luaDynamics)
    {
        L = luaL_newstate();
        luaL_openlibs(L); // might not be needed, performance penalty?
        lua_register(L, "usleep", usleep_lua);

        if (luaL_dofile(L, p->luaDynamicsFile))
            luaerror(L, "cannot run config. file: %s", lua_tostring(L, -1));
    }

    int luaquit = 0;

    size_t iter = 0;
    do{
        iter++;
        p->iter_final++;


        /* S -- Settings for this iteration
         *
         * Theses schemes typically includes some heuristics to avoid local minima
         */

        if(p->luaDynamics == 0)
        {
            fprintf(stderr, "luaDynamics == 0\n");
            exit(EXIT_FAILURE);
        }

        /* Settings from lua script */
        //     printf("Getting settings from lua script: %s\n", p->luaDynamicsFile);

        lua_getglobal(L, "getConfig"); /* function to be called */
        lua_pushnumber(L, iter); /* push arguments */
        lua_pushnumber(L, p->newx);

        /* do the call (2 arguments, 0 result) */
        if (lua_pcall(L, 2, 0, 0) != LUA_OK)
            luaerror(L, "error running function 'f': %s",
                     lua_tostring(L, -1));

        /* retrieve result */
        fconf.kDom = getglobfloat(L, "kDom");
        fconf.kVol = getglobfloat(L, "kVol");
        fconf.kInt = getglobfloat(L, "kInt");
        fconf.kRad = getglobfloat(L, "kRad");
        p->compress = getglobfloat(L, "kCom");
        Fb = getglobfloat(L, "fBrown");
        luaquit = getglobint(L, "quit");
        fconf.dInteraction = fconf.r0 * getglobfloat(L, "dInteraction");
        if(p->verbose > 10)
        {
            printf("kInteraction: %f\n", fconf.kInt);
            printf("dInteraction: %f\n", fconf.dInteraction);
        }

        /* Start of Molecular Dynamics
         * 2.1 Gradient of functional */
        grad3(X,
              p->N,
              p->R,
              p->I,
              g,
              &fconf);

        /* Cap gradient for stability */
        for(size_t kk = 0; kk<3*p->N; kk++)
        {
            if(fabs(g[kk]) > 1.0)
            {
                g[kk] = copysign(1.0, g[kk]);
            }
        }

        if(p->compress > 0)
        {
            comforce(p, X, g);
        }

        /* 2.2 Brownian force */
        if(Fb<0)
            Fb = 0;

        if(Fb>0)
        {
            for(size_t kk = 0; kk < p->N; kk++)
            {
                double d[] = {0,0,0};
                d[0] = normal();
                d[1] = normal();
                d[2] = normal();
#ifndef NDEBUG
                if(norm3(d)>10)
                {
                    printf("Unusual high norm of brownian: %f\n", norm3(d));
                    exit(-1);
                }
#endif
                //rand3d(&d[0]); // A unit length 3D direction
                //And give it a magnitude
                //double Frnd = 0.5*((double) rand() / (double) RAND_MAX - 0.5);
                for(size_t idx =0 ; idx<3; idx++)
                {
                    //              printf("%f -- ", norm3(g+3*kk));
                    g[3*kk+idx] += Fb*Frnd*d[idx];
                    //            printf("%f\n", norm3(g+3*kk));
                    // g[3*kk+idx] = (1-Fb)*g[3*kk+idx] + Fb*Frnd*d[idx];
                }
            }
        }

        /* 2.3 Dampening */
        /* Estimate velocities, note: per component */

        for(size_t pp = 0 ; pp < 3*p->N ; pp++)
        {
            v[pp] = (X[pp] - Xm[pp]) / (2.0 * dt);
        }

        for(size_t kk = 0; kk < 3*p->N; kk++)
        {
            g[kk] = g[kk] + damp*v[kk];
        }

        /* 3. Update X */
        for(size_t pp = 0 ; pp < 3*p->N ; pp++)
        {
            double xt = X[pp];
            X[pp] = 2.0*X[pp] - Xm[pp] - g[pp]*pow(dt,2);
            Xm[pp] = xt; // Update Xm to reflect the previous X-value
        }
        /* End of molecular dynamics */

        if(iter % 1500 == 0 || iter == maxiter-1)
        { showStep = 1; } else { showStep = 0; }

        if(showStep == 1)
        {
            /* Calculate gradient 2-norm */
            gnorm = 0;
            for(size_t kk = 0; kk < 3*p->N; kk++)
            {
                gnorm += pow(v[kk], 2);
            }
            gnorm = sqrt(gnorm);

            error = err3(X,
                         p->N,
                         p->R,
                         p->I,
                         &fconf);
            logwrite(p, 2, "    Iter: %6zu, E: %e, ||G||: %e\n",
                     iter, error, gnorm);
            fflush(p->logf);
        }

    } while( (iter < maxiter) && (run == 1) && (luaquit == 0));

    // TODO this is already calculated at at end of the last iteration.
    // also DRY.

    /* At final step, report back */
    double errorFinal = err3(X,
                             p->N,
                             p->R,
                             p->I,
                             &fconf);

    /* The gradient is a 3*p->N-dimensional vector, we return the 2-norm
       as the grad_final */
    double gradFinal = 0;
    for(size_t kk = 0; kk < 3*p->N; kk++)
    {
        gradFinal += pow(v[kk], 2);
    }

    gradFinal = sqrt(gradFinal);

    clock_gettime(CLOCK_MONOTONIC, &tb);

    p->err_final = errorFinal;
    p->grad_final = gradFinal;
    p->time_final = clockdiff(&ta, &tb);

    free(v);
    free(g);
    free(Xm);

    if(fconf.Es != NULL)
    {
        free(fconf.Es);
    }

    return 0;
}

/** @brief Report on success to log and screen
 *
 * @todo Time in ISO format.
 */
void param_summary(optparam * p, const double * restrict X)
{
    assert(p->N>0);

    logwrite(p, 1, "\n");
    logwrite(p, 1, " >> Optimization summary:\n");
    logwrite(p, 1, "    Final iterations: %zu\n", p->iter_final);
    logwrite(p, 1, "    Total time: %zu s\n", p->time_final);
    logwrite(p, 1, "    Final gradient norm: %e\n", p->grad_final);
    logwrite(p, 1, "    Final total error: %e\n", p->err_final);

    // X: mean, max, min
    double mex = 0, mey = 0, mez = 0;
    double md = 10e99;
    double mix = md, miy = md, miz = md;    /* min, x, y, z */
    double max = -md, may = -md, maz = -md; /* max, x, y, z */
    double mer = 0, mir = 10e99, mar = 0;   /* mean, min, max of radius */

    for(size_t kk = 0 ; kk<p->N; kk++)
    {
        mex += X[3*kk];
        mey += X[3*kk+1];
        mez += X[3*kk+2];
        mix = dmin(mix, X[3*kk]);
        miy = dmin(miy, X[3*kk+1]);
        miz = dmin(miz, X[3*kk+2]);
        max = dmax(max, X[3*kk]);
        may = dmax(may, X[3*kk+1]);
        maz = dmax(maz, X[3*kk+2]);
        double pr = norm3(X+3*kk);
        mar = dmax(mar, pr);
        mir = dmin(mir, pr);
        mer += pr;
    }
    mex /= p->N;
    mey /= p->N;
    mez /= p->N;
    mer /= p->N;

    logwrite(p, 2,  "\n");
    logwrite(p, 2,  ">> Structure summary:\n");
    logwrite(p, 2,  "              X       Y       Z       R\n");
    logwrite(p, 2,  "   Max:  % .3f, % .3f, % .3f, % .3f\n", max, may, maz, mar);
    logwrite(p, 2,  "   Mean: % .3f, % .3f, % .3f, % .3f\n", mex, mey, mez, mer);
    logwrite(p, 2,  "   Min:  % .3f, % .3f, % .3f, % .3f\n", mix, miy, miz, mir);

    if(run == 0)
    {
        logwrite(p, 0, "abnormal exit (Ctrl+c was pressed?)\n");
    }

}

/** @brief Read contact pairs from a binary file
    Sets p->I (the contacts) and p->NI (number of contact pairs)
    @return - Nothing, but aborts the program on failure.
*/
void param_read_contact_pairs(optparam * p)
{
    logwrite(p, 1, "Reading pairwise interactions from %s\n",
             p->contact_pairs_file);
    uint64_t nCP = 0;
    p->I = contact_pairs_read(p->contact_pairs_file, &nCP);
    if(p->I == NULL)
    {
        fprintf(stderr, "%s/%d Failed to read contact pairs from %s\n",
                __FILE__, __LINE__, p->contact_pairs_file);
        exit(EXIT_FAILURE);
    }
    p->NI = nCP;
    logwrite(p, 1, "Read %lu contacts pairs\n", nCP);
    return;
}

/** @brief Read binary indication matrix.

    Read the uint8 matrix, typically W.uint8, and set p->I to the
    contact pairs that can be found,
    p->I = [ a1, a2
    b1, b2,
    ...
    ]
    such that a1 < a2, b1 < b2 etc.
    also set p->NI
    also set p->N
    depreciated. Use param_read_contact_pairs
*/
void param_readW(optparam * p)
{

    logwrite(p, 1, "Reading pairwise interactions from %s\n", p->wfname);

    size_t fsize = 0;
    uint8_t * W = wio_read(p->wfname, &fsize);
    if(fsize == 0)
    {
        printf("%s/%d Failed reading %s!\n",
               __FILE__,  __LINE__, p->wfname);
        exit(EXIT_FAILURE);
    }

    logwrite(p, 1, " W = [ %d, %d, %d, %d, ..., %d, %d, %d, %d]\n",
             (int) W[0],
             (int) W[1],
             (int) W[2],
             (int) W[3],
             W[fsize-4],
             W[fsize-3],
             W[fsize-2],
             W[fsize-1]);

    // Set the number of elements
    p->N = (size_t) sqrt(fsize);
    assert(p->N*p->N == fsize);
    logwrite(p, 1, " %zu points\n", p->N);

    // Initiate the list with pairs in contact
    // I.e. p->I and p->NI

    // Count them
    size_t nPairs = 0;
    for(size_t kk = 0; kk<p->N; kk++)
    {
        for(size_t ll = kk+1; ll<p->N; ll++)
        {
            if(W[kk+ll*p->N] == 1)
            {
                nPairs++;
            }
        }
    }

    logwrite(p, 1, " Found %zu interaction pairs\n", nPairs);

    if(nPairs > 12*p->N)
    {
        logwrite(p, 0, "more than 12 pairs per bead in average!\n");
    }

    p->NI = nPairs;
    p->I = malloc(nPairs*2*sizeof(uint32_t));
    assert(p->I != NULL);

    // Construct the list
    size_t writepos = 0;
    for(size_t kk = 0; kk<p->N; kk++)
        for(size_t ll = kk+1; ll<p->N; ll++)
        {
            if(W[kk+ll*p->N] == 1)
            {
                p->I[writepos++] = kk;
                p->I[writepos++] = ll;
            }
        }
    printf(" ok\n");
    free(W);
    return;
}

double * init_X(optparam * p)
{
    double * X = NULL;

    p->newx = 0;
    if(p->xfname != NULL)
    {
        X = malloc(3*p->N*sizeof(double));
        assert(X != NULL);
        if(param_readX(p, X) != 0)
        {
            logwrite(p, 2, "Could not open x-file, starting from random\n");
            free(X);
            X=NULL;
        }
    }

    if(X==NULL)
    {
        p->newx = 1;
        logwrite(p, 1, "Using random initialization for X\n");
        srand(p->rseed);
        X = malloc(3*p->N*sizeof(double));
        assert(X != NULL);
        elli * E = p->E;
        if(E == NULL){
            E = elli_new(1, 1, 1);
        }

        for(size_t kk = 0; kk< p->N; kk++)
        {
            int accepted = 0;
            while(accepted == 0)
            {
                for(int idx =0; idx<3; idx++)
                {
                    X[3*kk+idx] = 2.0*(rand()/(double) RAND_MAX-.5);
                }
                if(elli_getScale(E, X+3*kk) <= 1)
                { accepted = 1;}
            }
        }

        logwrite(p, 1, "X[0] = %f\n", X[0]);
    }
    return X;
}

int param_readL(optparam * p)
{
    /* Read label matrix pointed to by p->lfname
     */

    if(p->lfname == NULL)
    {
        logwrite(p, 1, "No L-file specified\n");
        return EXIT_FAILURE;
    }

    logwrite(p, 1, "Reading L-labels from %s\n", p->lfname);

    // Try to read as binary
    FILE * f = fopen(p->lfname, "rb");
    if(f == NULL)
    {
        fprintf(stderr, "Unable to read %s\n", p->lfname);
        exit(EXIT_FAILURE);
    }
    fseek(f, 0, SEEK_END); // seek to end of file
    size_t fsize = ftell(f); // get current file pointer
    fseek(f, 0, SEEK_SET); // seek back to beginning of file

    logwrite(p, 1, "As uint8_t, %s %zu numbers (%zu bytes)\n", p->lfname,
             fsize/sizeof(uint8_t), fsize);

    p->N = fsize/sizeof(uint8_t);
    p->L = malloc(p->N*sizeof(double));
    assert(p->L != NULL);
    size_t nread = fread(p->L, sizeof(uint8_t), p->N, f);
    if(nread != p->N)
    {
        fprintf(stderr, "Unable to read from %s\n", p->lfname);
        exit(EXIT_FAILURE);
    }
    fclose(f);

    logwrite(p, 2, "L = [%u, %u, ..., %u]\n", p->L[0], p->L[1], p->L[p->N-1]);

    return 0;
}

int param_readR(optparam * p)
{
    /* Read GPSeq radius values as binary double.
     * If that does not work, try as text, one value per line
     */
    if(p->rfname == NULL)
    {
        logwrite(p, 1, "No radial preferences to read\n");
        return EXIT_FAILURE;
    }

    logwrite(p, 1, "Reading R-values from %s\n", p->rfname);
    size_t nbytes = 0;
    p->R = (double *) wio_read(p->rfname, &nbytes);
    printf("Read %zu doubles from %s\n", nbytes/sizeof(double), p->rfname);

    if(2*nbytes/sizeof(double) == p->N)
    {
        printf("Found half as many values as beads, assuming diploid and duplicating data\n");
        if(p->N%2 == 1)
        {
            printf("ERROR: N%%2 == 1\n");
            exit(-1);
        }

        p->R = realloc(p->R, p->N*sizeof(double));
        if(p->R == NULL)
        {
            printf("ERROR: out of memory\n");
            exit(-1);
        }
        for(size_t pp = 0; pp < p->N/2; pp++)
        {
            p->R[pp + p->N/2] = p->R[pp];
        }
        //memcpy(p->R+p->N/2*sizeof(double), p->R, p->N/2*sizeof(double));
    }


    if(!( (nbytes/sizeof(double) == p->N) || (nbytes/sizeof(double) == p->N/2) ))
    {
        printf("Error: Can't make sense of %s, expected %zu or %zu bytes but got %zu\n", p->rfname, 4*p->N, 2*p->N, nbytes);
        exit(-1);
    }

    size_t nInf = 0;
    for(size_t kk = 0; kk<p->N; kk++)
    {
        if(!isfinite(p->R[kk]))
        {
            nInf++;
        }
    }

    printf("%s contains %zu non-finite values which will be ignored.\n", p->rfname, nInf);

    return EXIT_SUCCESS;
}

/* Write coordinates to disk, also write the column names
   to the log file  */
void param_dumpX(optparam * p, double * X)
{
    size_t N = p->N;

    int write_R = 0;
    if(p->R != NULL)
        write_R = 1;

    if(run == 1)
    {
        if(p->verbose > 1)
        {
            logwrite(p, 1, "Writing final structure to %s\n", p->xoutfname);
        }
    }
    else {
        p->xoutfname = realloc(p->xoutfname, sizeof(char)*(strlen(p->xoutfname)+5));
        strcat(p->xoutfname, ".0");
        //    sprintf(p->xoutfname, "%s.0", p->xoutfname);

        logwrite(p, 1, "Writing non-finished structure to: %s\n", p->xoutfname);
    }

    if(p->verbose > 1)
    {
        logwrite(p, 1, "Columns: x, y, z, r");

        if(write_R)
        {
            logwrite(p, 1, ", R");
        }
        logwrite(p, 1, "\n");
    }

    FILE * f = fopen(p->xoutfname, "w");

    for(size_t kk = 0; kk<N; kk++)
    {
        double radius;
        if(p->E == NULL)
        {
            radius = norm3(X+3*kk);
        } else {
            radius = elli_getScale(p->E, X+3*kk);
        }

        fprintf(f, "%f, %f, %f, %f", X[3*kk], X[3*kk+1], X[3*kk+2], radius);
        if(write_R)
        {
            fprintf(f, ", %f", p->R[kk]);
        }

        fprintf(f, "\n");
    }
    fclose(f);
    return;
}

void param_show(optparam * p, FILE * f)
{
    double volocc = p->N*4.0/3.0*M_PI*pow(p->r0,3) / (4.0/3.0*M_PI);

    if(p->E != NULL)
    {
        volocc = p->N*4.0/3.0*M_PI*pow(p->r0,3) / elli_vol(p->E);
    }


    fprintf(f, "\n");
    fprintf(f, " >> Parameters:\n");
    fprintf(f, "    Problem size: %zu points (%zu variables)\n", p->N, 3*p->N);
    fprintf(f, "    - Model parameters\n");
    if(p->luaDynamics == 0)
    {
        fprintf(f, "    Spring constant volume exclusion (repulsion): %f\n", p->kVol);
        fprintf(f, "    Spring constant bead interactions (attraction): %f\n", p->kInt);
        fprintf(f, "    Spring constant domain confinement: %f\n", p->kDom);
        fprintf(f, "    Spring constant radial positioning: %f\n", p->kRad);
        fprintf(f, "    Max iterations: %zu\n", p->maxiter);
        fprintf(f, "    Max time: %zu (s)\n", p->maxtime);
        if(p->compress == 1)
        {
            fprintf(f, "    Compression: on\n");
        } else {
            fprintf(f, "    Compression: off\n");
        }
    }
    if(p->luaDynamics == 1)
    {
        fprintf(f, "    Dynamics program: %s\n", p->luaDynamicsFile);
    }
    fprintf(f, "    Bead radius %f (vol. occ. %f)\n", p->r0, volocc);
    fprintf(f, "    - Optimization parameters\n");
    //fprintf(f, "# Optimization: simulated annealing / molecular dynamics\n");
    fprintf(f, "    random seed: %zu\n", p->rseed);
    fprintf(f, "    write compresseed cmm: %d\n", p->cmmz);

    if(p->wfname == NULL)
    {
        fprintf(f, "    W file: -not specified\n");
    } else
    {
        fprintf(f, "    W file: %s\n", p->wfname);
    }

    if(p->contact_pairs_file == NULL)
    {
        fprintf(f, "   Contact Pairs file not specified (REQUIRED!)\n");
    } else {
        fprintf(f, "   Contacts Pairs File: %s\n", p->contact_pairs_file);
    }

    if(p->xfname == NULL)
    {
        fprintf(f, "    X file: -not specified-\n");
    }
    else
    {
        fprintf(f, "    X file: %s\n", p->xfname);
    }

    if(p->lfname == NULL)
    {
        fprintf(f, "    L file: -not specified-\n");
    }
    else
    {
        fprintf(f, "    L file: %s\n", p->lfname);
    }

    if(p->rfname == NULL)
    {
        fprintf(f, "    R file: -not specified-\n");
    }
    else
    {
        fprintf(f, "    R file: %s\n", p->rfname);
    }

    if(p->ofoldername == NULL)
    {
        fprintf(f, "    output folder not specified\n");
    }
    else
    {
        fprintf(f, "    output folder: %s\n", p->ofoldername);
    }

    fprintf(f, "    verbose level: %d\n", p->verbose);
    fprintf(f, "\n");
}


/* do-nothing callback */
static int usleep_lua_NULL(lua_State * L)
{
    /* get number of arguments */
    if(lua_gettop(L) != 1)
    {
        lua_pushstring(L, "Incorrect number of arguments to 'usleep'");
        lua_error(L);
        return 0;
    }

    if (!lua_isnumber(L, 1))
    {
        lua_pushstring(L, "Incorrect argument to 'usleep'");
        lua_error(L);
        return 0;
    }

    /* return the number of results */
    return 0;
}

/* Write out the behavior of the lua script as a table */
void dump_lua_dynamics(char * luafile)
{
    printf("newx, iter, kDom, kVol, kInt, dInt_rel, kRad, compress, Fb\n");
    for(int newx = 1; newx >= 0; newx--)
    {
        lua_State * L = luaL_newstate();
        luaL_openlibs(L);
        /* Do nothing when usleep is called */
        lua_register(L, "usleep", usleep_lua_NULL);

        if (luaL_dofile(L, luafile))
            luaerror(L, "cannot run config. file: %s", lua_tostring(L, -1));

        size_t iter = 0;
        int luaquit = 0;
        char luafun[] = "getConfig";
        do{
            iter++;
            lua_getglobal(L, luafun); /* function to be called */
            lua_pushnumber(L, iter); /* push arguments */
            lua_pushnumber(L, newx);
            /* do the call (2 arguments, 0 result) */
            if (lua_pcall(L, 2, 0, 0) != LUA_OK)
                luaerror(L, "error running function '%s': %s",
                         luafun,
                         lua_tostring(L, -1));

            double kDom = getglobfloat(L, "kDom");
            double kVol = getglobfloat(L, "kVol");
            double dInt = getglobfloat(L, "dInteraction");
            double kInt = getglobfloat(L, "kInt");
            double kRad = getglobfloat(L, "kRad");
            double compress = getglobfloat(L, "kCom");
            double Fb = getglobfloat(L, "fBrown");
            luaquit = getglobint(L, "quit");
            //fconf.dInteraction = fconf.r0 * getglobfloat(L, "dInteraction");
            printf("%d, %zu, %f, %f, %f, %f, %f, %f, %f\n",
                   newx, iter, kDom, kVol, kInt, dInt, kRad, compress, Fb);
        } while(luaquit == 0);
        lua_close(L);
    }
    return;
}

void usage()
{
    printf("Note that mflock typically is run via chromflock\n"
           "see man chromflock for more details\n");
    printf("1. Required arguments:\n");
    printf(" --contact-pairs <file>\n\tFile with contact pairs (uint32).\n");
    printf(" -L <file>, --lFile <file>\b\tlabel matrix (uint8_t)\n");
    printf(" --dconf mflock.lua\n\tSpecify lua script to use for dynamics\n");
    printf("2. Optional arguments:\n");
    printf(" -x <file>, --xFile <file>\n\tfile with initial coordinates\n");
    printf(" -s N, --seed N\n\tseed for random number generator (defaults to random)\n");
    printf(" -r <file>, --rFile <file>\n\tfile with wanted radii (in conjunction with -G)\n"
           "\tall files should be encoded as 64-bit floats\n");
    printf(" -n N, --maxiter N\n\tmaximum number of iterations\n");
    printf(" -t N, --maxtime N\n\t time budget in seconds\n");
    printf("2.1 Forces\n");
    printf(" --dconf-show mflock.lua\n\tPrint the dynamics as a table and quit.\n");
    printf("\tNo additional arguments required or used\n");
    printf("2.2 Geometry\n");
    printf(" -R r0, --radius r0\n\tbead radius (sphere has radius 1)\n");
    printf(" -Q vq, --volq vq\n\tvolume quotient (beads/domain)\n"
           "    For ellipsoidal geometry, set 1=ea>=eb>=ec>0\n");
    printf(" -A ea, --ea ea\n\tMajor axis of ellipsoid\n");
    printf(" -B eb, --eb eb\n\tSecond axis of ellipsoid\n");
    printf(" -C ec, --ec ec\n\tThird axis of ellipsoid\n");
    printf("2.3 General\n");
    printf(" -v level, --verbose level\n\tverbosity lowest=0\n");
    printf(" -o name, --outFolder name\n\twhere to store results\n");
    printf(" -z, --cmmz\n\twrite compressed cmm files (gzip)\n");
    printf(" -c, --compress\n\tenable chromosome compression (only for -D)\n");
    printf(" -a, --live\n\tenable live monitoring (only when compiled with SDL)\n");
    printf(" -h, --help\n\tshow this help message. For more info see 'man mflock'\n");
    printf(" -d, --defaults\n\tshow default settings for the parameters\n");
    printf("\n");
    return;
}

#define ARGS_OK 0
#define ARGS_ERR 1
#define ARGS_QUIT 2

int argparsing(optparam * p, int argc, char ** argv)
{

    /* Specifications of ellipsoid */
    double ea = -1;
    double eb = -1;
    double ec = -1;

    struct option longopts[] = {
        { "version",      no_argument,       NULL,   'i' },
        { "help",         no_argument,       NULL,   'h' },
        // Data
        { "wFile",        required_argument, NULL,   'w' },
        { "contact-pairs", required_argument, NULL,  'p' },
        { "xFile",        required_argument, NULL,   'x' },
        { "rFile",        required_argument, NULL,   'r' },
        { "labels",        required_argument, NULL,   'L' },
        { "outFolder",    required_argument, NULL,   'o' },
        // Settings
        { "maxiter",      required_argument, NULL,   'n' },
        { "maxtime",      required_argument, NULL,   't' },
        { "seed",         required_argument, NULL,   's' },
        { "verbose",      required_argument, NULL,   'v' },
        { "live",         no_argument, NULL,   'a' },
        { "cmmz",         no_argument, NULL,   'z' },
        // Settings / Forces
        { "kVol",        required_argument, NULL,   'V' },
        { "kDom",        required_argument, NULL,   'D' },
        { "kInt",        required_argument, NULL,   'I' },
        { "kRad",        required_argument, NULL,   'G' },
        // Settings / Geometry
        { "radius",        required_argument, NULL,   'R' },
        { "vq",          required_argument, NULL,   'Q' },
        { "ea",            required_argument, NULL,   'A' },
        { "eb",            required_argument, NULL,   'B' },
        { "ec",            required_argument, NULL,   'C' },
        // Settings / Optimization
        { "compress",       no_argument,       NULL,   'c' },
        { "defaults",       no_argument,       NULL,   'd' },
        // Lua program for dynamics
        { "dconf",         required_argument,    NULL,   'l' },
        { "dconf-show",    required_argument,    NULL,   'M' },
        { NULL,           0,                 NULL,   0   }
    };

    int ch;
    while((ch = getopt_long(argc, argv,
                            "A:B:C:w:x:r:n:t:V:I:S:R:G:v:o:p:hMs:L:zcadQ:l:",
                            longopts, NULL)) != -1)
    {
        switch(ch) {
        case 'i':
            printf("mflock (chromflock version %s)\n", cf_version);
            printf("Build date: %s, %s\n", __DATE__, __TIME__);
            printf("GIT HASH: %s\n", GIT_VERSION);
            printf("Compiler: %s\n", CC_VERSION);
            return ARGS_QUIT;
        case 'd':
            printf("Defaults:\n");
            param_show(p, stdout);
            return ARGS_QUIT;
        case 'w':
            if(p->wfname != NULL)
            {free(p->wfname);}
            p->wfname = malloc(strlen(optarg)+1);
            assert(p->wfname != NULL);
            strcpy(p->wfname, optarg);
            break;
        case 'x':
            p->xfname = malloc(strlen(optarg)+1);
            assert(p->xfname != NULL);
            strcpy(p->xfname, optarg);
            break;
        case 'r':
            p->rfname = malloc(strlen(optarg)+1);
            assert(p->rfname != NULL);
            strcpy(p->rfname, optarg);
            break;
        case 'L':
            p->lfname = malloc(strlen(optarg)+1);
            assert(p->lfname != NULL);
            strcpy(p->lfname, optarg);
            break;
        case 'n':
            p->maxiter = atol(optarg);
            break;
        case 'p':
            free(p->contact_pairs_file);
            p->contact_pairs_file = strdup(optarg);
            break;
        case 't':
            p->maxtime = atol(optarg);
            break;
        case 'V':
            p->kVol = atof(optarg);
            break;
        case 'S':
            p->kDom = atol(optarg);
            break;
        case 'I':
            p->kInt = atof(optarg);
            break;
        case 'G':
            p->kRad = atof(optarg);
            break;
        case 'R':
            p->r0 = atof(optarg);
            break;
        case 'Q':
            p->volq = atof(optarg);
            break;
        case 's':
            p->rseed = atol(optarg);
            break;
        case 'v':
            p->verbose = atoi(optarg);
            break;
        case 'o':
            free(p->ofoldername);
            p->ofoldername = malloc(strlen(optarg)+2);
            assert(p->ofoldername != NULL);
            strcpy(p->ofoldername, optarg);
            break;
        case 'l':
            p->luaDynamics = 1;
            p->luaDynamicsFile = malloc(strlen(optarg) + 1);
            assert(p->luaDynamicsFile != NULL);
            strcpy(p->luaDynamicsFile, optarg);
            break;
        case 'M':
            dump_lua_dynamics(optarg);
            return ARGS_QUIT;
        case 'h':
            return(1);
        case 'z':
            p->cmmz = 1;
            break;
        case 'c':
            p->compress = 1;
            break;
        case 'a':
            p->liveView = 1;
            break;
        case 'A':
            ea = atof(optarg);
            break;
        case 'B':
            eb = atof(optarg);
            break;
        case 'C':
            ec = atof(optarg);
            break;
        default:
            return ARGS_ERR;
        }
    }

    if(p->luaDynamics == 0)
    {
        fprintf(stderr, "--dconf not specified\n");
        exit(EXIT_FAILURE);
    }

    /* Make sure that the outfolder ends with a path_separator.
       Note: it should already be allocated to have room for it if missing */
    if(p->ofoldername)
    {
#ifdef _WIN32
        fprintf(stderr, "TODO: Non-portable section %s %s\n",
                __FILE__, __LINE__);
        exit(EXIT_FAILURE);
#endif
        if(p->ofoldername[strlen(p->ofoldername)] != '/' )
        {
            p->ofoldername[strlen(p->ofoldername)+1] = '\0';
            p->ofoldername[strlen(p->ofoldername)] = '/';
        }
    }

    if(p->contact_pairs_file == NULL)
    {
        fprintf(stderr, "WARNING: Contacts file not set\n");

        if(p->wfname == NULL)
        {
            printf("ERROR: 'W' file not set (-w) \n");
            if(0){
                printf("Generate one in MATLAB with:\n");
                printf("A = zeros(150,150)\n");
                printf("A(2:size(A,1)+1:end) = 1;\n");
                printf("A = A + A';\n");
                printf("A8 = uint8(A);\n");
                printf("fout = fopen('A150.dat', 'wb');\n");
                printf("fwrite(fout, A, 'uint8');\n");
                printf("fclose(fout)\n");
                printf("\n\n");
            }
            return ARGS_ERR;
        }
    }

    int efail = 0;
    if(ea > 0 || eb > 0 || ec > 0)
    {
        if(ea>= eb && eb >= ec)
        {
            if(ea == 1 && ec > 0)
            {
                //        printf("Using ellipsoidal geometry!\n");
                p->E = elli_new(ea, eb, ec);
                // elli_show(p->E);
            } else {efail = 1;}
        } else {efail = 1;}
    }
    if(efail)
    {
        printf("ERROR: Ellipsoidal geometry requires 1=ea>=eb>=ec>0\n\n");
        return ARGS_ERR;
    }

    return ARGS_OK;
}

optparam *  param_alloc(void)
{
    optparam * p = calloc(1, sizeof(optparam));
    assert(p != NULL);

    struct timespec ts;
    //  timespec_get(&ts, TIME_UTC);
    clock_gettime(CLOCK_REALTIME, &ts);

    p->rseed = time(NULL)*getpid()*ts.tv_nsec;

    p->maxiter = 1000000; // iterations
    p->maxtime = 60*60*10; // seconds

    p->grad_final = 0;
    p->iter_final = 0;
    p->time_final = 0;

    p->kVol = 1;
    p->kInt = 1;
    p->kDom = 1;
    p->kRad = 0;
    p->r0 = -1;
    p->volq = 0.2;

    p->N = 0;
    p->I = NULL;
    p->NI = 0;
    p->R = NULL;
    p->verbose = 1;
    p->compress = 0;
    p->cmmz = 0;
    p->liveView = 0;

    p->L = NULL;

    p->E = NULL;
    p->newx = 1;

    // Create a suggestion for the output folder
    if(p->ofoldername == NULL)
    {

        p->ofoldername = malloc(1024*sizeof(char));
        assert(p->ofoldername != NULL);

        size_t fn = 1;
        int okfolder = 0;

        struct stat st;
        memset(&st, 0, sizeof(struct stat));

        while(okfolder == 0)
        {
            sprintf(p->ofoldername, "sf_%05zu/", fn);
            fn++;
            if(stat(p->ofoldername, &st) == -1)
                okfolder = 1;
        }
    }

    fflush(stdout);
    return p;
}

int param_init(optparam * p, int argc, char ** argv)
{

    struct stat st;
    memset(&st, 0, sizeof(struct stat));

    if(stat(p->ofoldername, &st) == -1)
    {
        if(p->verbose >= 1){
            printf("Creating output folder: %s\n", p->ofoldername);
        }
        int folderok = mkdir(p->ofoldername, 0770);
        if(folderok != 0)
        {
            printf("Could not create output folder\n");
            exit(-1);
        }

    } else {
        if(p->verbose >= 2){
            printf("Output folder did already exist.\n");
        }
    }

    char * fullofolder = realpath(p->ofoldername, NULL);
    if(fullofolder == NULL)
    {
        printf("Error: The output folder '%s' can't be understood by realpath\n", p->ofoldername);
        exit(-1);
    }

    // Create output folder
    if(p->verbose >= 1) {
        printf("Output folder: %s\n", fullofolder);
    }
    free(fullofolder);


    // Set names of output files
    p->xoutfname = malloc(1024*sizeof(char));
    assert(p->xoutfname != NULL);
    sprintf(p->xoutfname, "%s%s", p->ofoldername, "coords.csv");

    p->logfname = malloc(1024*sizeof(char));
    assert(p->logfname != NULL);
    sprintf(p->logfname, "%s%s", p->ofoldername, "log.txt");

    // Open log
    p->logf = fopen(p->logfname, "a");
    assert(p->logf != NULL);

    if(p->logf== NULL)
    {
        printf("Failed to open log file\n");
        exit(1);
    }

    char * time_str = cf_timestr();
    fprintf(p->logf, "\nmflock started: %s\n", time_str);
    free(time_str);

    fprintf(p->logf, "CMD: ");
    for(int kk = 0; kk<argc; kk++)
    {
        fprintf(p->logf, "%s ", argv[kk]);
    }
    fprintf(p->logf, "\n");
    fflush(p->logf);

    return 0;
}


void param_free(optparam * p)
{
    free(p->R);
    free(p->I);
    free(p->L);
    free(p->wfname);
    free(p->contact_pairs_file);
    free(p->lfname);
    free(p->rfname);
    free(p->xfname);
    free(p->xoutfname);
    free(p->ofoldername);
    free(p->luaDynamicsFile);
    free(p->logfname);
    free(p->E);

    free(p);
    return;
}


void param_validate(optparam * p)
{
    if(p->N == 0)
    {
        fprintf(stderr, "p->N = 0\nCheck the size of the W-matrix\n");
        exit(1);
    }

    if(p->logf == NULL)
    {
        fprintf(stderr, "p->logf = NULL\nThis indicates that the log file could not be opened\n");
        exit(1);
    }

    if(p->r0 < 0)
    {
        double vq = p->volq;
        if(p->verbose > 1)
        {
            printf("Bead radius not set, setting total bead volume to %f%%\n", 100*vq);
        }
        double Vd = 4.0/3.0*M_PI;
        if(p->E != NULL)
        {
            Vd = elli_vol(p->E);
        }
        p->r0 = cbrt( 3.0*vq*Vd / (4.0*p->N*M_PI) );

#ifndef NDEBUG
        double Vbeads = p->N*pow(p->r0,3)*M_PI*4.0/3.0;
        printf("Vd: %f, Vbeads: %f, Vbeads/Vd = %f\n", Vd, Vbeads, Vbeads/Vd);
        assert(fabs(Vbeads/Vd - vq)<1e-5);
#endif
    }
}

void logwrite(optparam * p, int level, const char *fmt, ...)
{

    if(p->verbose >= level)
    {
        va_list args, args2;
        va_start(args, fmt);
        va_copy(args2, args);
        fprintf(stdout, "    ");
        if(level == 0) {
            fprintf(stdout, "WARNING: "); }
        vfprintf(stdout, fmt, args);
        va_end(args);
        if(p->logf != NULL)
        {
            vfprintf(p->logf, fmt, args2);
        } else {
            printf("Log file not available\n");
        }
        va_end(args2);
    }
    return;
}


int param_readX(optparam * p, double * X)
{

    //  fprintf(stdout, "Reading X-data from %s\n", p->xfname);
    FILE * f = fopen(p->xfname, "r");
    if(f == NULL)
    {
        printf("Can't open %s\n", p->xfname);
        return -1;
    }

    char * line = malloc(1024*sizeof(char));
    assert(line != NULL);
    size_t len = 0;

    char delim[] = ",";
    for(size_t ll = 0; ll<p->N; ll++)
    {
        int read = getline(&line, &len, f);
        if(read == -1)
        {
            printf("Failed to read line %zu\n", ll+1);
            return -1;
        }
        char *ptr = strtok(line, delim);
        X[3*ll] = atof(ptr);
        ptr = strtok(NULL, delim);
        X[3*ll+1] = atof(ptr);
        ptr = strtok(NULL, delim);
        X[3*ll+2] = atof(ptr);
    }

    free(line);
    return 0;
}


void optimize(double * restrict X, optparam * restrict p)
{
    if(p->luaDynamics == 0)
    {
        double Fb = .7;
        while(Fb>0.1)
        {
            dynamic(X, p, Fb);
            Fb = Fb*.7;
        }
        dynamic(X, p, 0);

        if(p->compress == 1)
        {
            printf("Running again with compression off (relaxation)\n");
            p->compress = 0;
            dynamic(X, p, 0);
        }
    } else {
        dynamic(X, p, 0); /* Default way */
    }

}


void * solve_t(void * args)
{
    solve_t_args * sargs = (solve_t_args *) args;
    optimize(sargs->X, sargs->p);
    return NULL;
}


int main(int argc, char ** argv)
{
    time_t starttime, nowtime;
    time(&starttime);

    optparam * p = param_alloc(); // Set default options

    switch(argparsing(p, argc, argv))
    {
    case ARGS_OK:
        break;
    case ARGS_ERR:
        usage();
        param_free(p);
        return EXIT_FAILURE;
        break;
    case ARGS_QUIT:
        param_free(p);
        return EXIT_SUCCESS;
        break;
    }

    param_init(p, argc, argv);

    if(p->contact_pairs_file != NULL)
    {
        param_read_contact_pairs(p);
    } else {
        if(p->wfname != NULL)
        {
            /* Read binary contact indication matrix */
            param_readW(p);
        } else {
            fprintf(stderr, "No contact pairs specified\n");
            exit(EXIT_FAILURE);
        }
    }

    /* Read radial preferences, if rfname is set */
    param_readR(p);

    /* Read chromosome labels */
    param_readL(p);

    /* Read coordinates from previous iteration, or generate
       random coordinates if not previous iteration is available
    */
    double * X = init_X(p);

    /* Double check that the parameters make sense and that
       everything necessary could be loaded. */
    param_validate(p);

    if(p->verbose>0)
    {
        param_show(p, stdout);
    }
    param_show(p, p->logf);

    /* Set up capturing of Ctrl+C for graceful exit */
    struct sigaction act;
    memset (&act, '\0', sizeof(act));
    act.sa_handler = stoprun;
    sigaction(SIGINT, &act, NULL);

    /* Start molecular dynamics */
    logwrite(p, 1, " >> Solving ... \n");

#ifdef SDL
    if(p->liveView == 1)
    {
        int quit = 0;

        pthread_t th;
        solve_t_args solve_arg;
        solve_arg.X = X;
        solve_arg.p = p;

        pthread_create(&th, // thread
                       NULL, // pthread_attrib_t
                       solve_t, // function
                       &solve_arg); // arg

        elli * E;
        if(p->E == NULL)
        {
            E = elli_new(1,1,1);
        } else {
            E = p->E;
        }

        liveview(X, p->L, p->N, &quit, p->r0, E);
        pthread_join(th, NULL);
    } else {
        optimize(X,p);
    }
#else
    optimize(X, p);
#endif
    //radpos(X, p, stdout);

    time(&nowtime);
    p->time_final = difftime(nowtime, starttime);

    /* Write a summary of how the simulations proceeded to the log file */
    param_summary(p, X);

    /* Write chimera file (for simple visualization) */
    if(p->cmmz == 1)
    {
        char * cmmfile = malloc(1024*sizeof(char));
        assert(cmmfile != NULL);
        sprintf(cmmfile, "%s/cmmdump.cmm.gz", p->ofoldername);

        cmmwritez(cmmfile, X, p->N, p->r0, p->I, p->NI, p->L);
        free(cmmfile);
    } else {
        char * cmmfile = malloc(1024*sizeof(char));
        assert(cmmfile != NULL);
        sprintf(cmmfile, "%s/cmmdump.cmm", p->ofoldername);
        cmmwrite(cmmfile, X, p->N, p->r0, p->I, p->NI, p->L);
        free(cmmfile);
    }

    /* Write coordinates to disk */
    // TODO: rename mflock_write_coords(mflock * m, double * X)
    param_dumpX(p, X);

    /* Write radial profile to disk ? */
    // Just wastes space, can as well be calculated on the fly.
    if(0) {
        char * rpfname = malloc(1024*sizeof(char));
        assert(rpfname != NULL);
        sprintf(rpfname, "%s/rad.csv", p->ofoldername);
        FILE * rpfile = fopen(rpfname, "w");
        assert(rpfile != NULL);
        radpos(X, p, rpfile);
        fclose(rpfile);
        free(rpfname);
    }

    char * time_str = cf_timestr();
    fprintf(p->logf, "\nmflock finished: %s\n", time_str);
    free(time_str);

    /* Close log file */
    fclose(p->logf);
    /* Free most data */
    param_free(p);
    /* Free coordinates */
    free(X);

    return EXIT_SUCCESS;
}
