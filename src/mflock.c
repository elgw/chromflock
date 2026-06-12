/**
 * @file mflock.c
 * @author Erik Wernersson
 * @date 2020-2023
 */

#include "mflock.h"
#include "mflock_private.h"

static bpos *  mflock_load_bead_apos(const char * fname, int * nbpos)
{
    if(npy_extension(fname))
    {
        int nconstraint = 0;
        bpos * apos =  load_bead_apos_from_npy(fname,
                                               &nconstraint);
        if(apos == NULL)
        {
            fprintf(stderr, "Failed to load absolute bead positions from %s\n", fname);
            exit(EXIT_FAILURE);
        }
        *nbpos = nconstraint;
        return apos;
    } else {
        fprintf(stderr,
                "Can only load absolute bead positions from npy files\n");
        exit(EXIT_FAILURE);
    }
    return NULL;
}

static double norm3d(const double * restrict X)
{
    double n = 0;
    for(size_t kk = 0; kk<3; kk++)
        n+=pow(X[kk], 2);
    return sqrt(n);
}


/* Used for communication with the optional visualization routines */
static volatile int run = 1;

/* show lua error string */
static void luaerror(lua_State *L, const char *fmt, ...)
{
    va_list argp;
    va_start(argp, fmt);
    vfprintf(stderr, fmt, argp);
    printf("\n");
    va_end(argp);
    lua_close(L);
    exit(EXIT_FAILURE);
}

/* get global variable by name from lua */
static double lua_get_float (lua_State *L,
                             const char *var)
{
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

/* get an integer from lua */
static int lua_get_int (lua_State *L, const char *var) {
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



static double dmax(double a, double b)
{
    if(a>b)
        return a;
    return b;
}

static double dmin(double a, double b)
{
    if(a<b)
        return a;
    return b;
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


/* Euclidean distance between two 3D-vectors - squared */
static double eudist3p2(const double * A, const double * B)
{
    return pow(A[0]-B[0], 2) + pow(A[1]-B[1], 2) + pow(A[2]-B[2], 2);
}


static void comforce(mflock_t * restrict p,
                     double * restrict G)
{
    const double * restrict X = p->beads;
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
    for(size_t kk = 0 ; kk < p->n_beads ; kk++)
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
    for(size_t kk = 0; kk < p->n_beads; kk++)
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


/* Here is the main loop */
static int
mflock_dynamics(mflock_t * restrict p)
{
    double * X = p->beads;
    /* set up fast prng in normal.h before calling normal() */
    normal_setup();

    /* Prepare settings */
    size_t maxiter = p->maxiter;
    mflock_func_t fconf = {0};
    fconf.r0 = p->r0;
    fconf.E = NULL;
    fconf.Es = NULL;
    fconf.diploid = p->diploid;
    fconf.top_plane = 2.0;
    fconf.bottom_plane = -2.0;
    fconf.geometry = p->geometry;
    if(p->E != NULL)
    {
        /* add ellipse parameters otherwise sphere domain */
        fconf.E = p->E;
        fconf.Es = elli_new(
                            p->E->a - fconf.r0,
                            p->E->b - fconf.r0,
                            p->E->c - fconf.r0);
        fconf.geometry = MFLOCK_ELLIPSOID;
    }

    fconf.nIPairs = p->n_pairs; /* Only for err2 */

    /* Initialization */
    struct timespec ta, tb;
    clock_gettime(CLOCK_MONOTONIC, &ta);

    /* Gradient */
    double * restrict g = calloc(p->n_beads*3, sizeof(double));

    /* Velocity */
    double * restrict v = calloc(p->n_beads*3, sizeof(double));

    /* Set previous position, initially to be the current position,
       i.e. initial velocity will be 0 */
    double * restrict Xm = malloc(3*p->n_beads*sizeof(double));
    memcpy(Xm, X, 3*p->n_beads*sizeof(double));

    double error = 9e99;
    double gnorm = 9e99;
    double dt = 0.15;
    double damp = 0.5; /* dampening (.8) faster convergence than .4 */
    int showStep = 0;

    lua_State *L = luaL_newstate();
    luaL_openlibs(L);
    lua_register(L, "usleep", usleep_lua);

    if (luaL_dofile(L, p->luaDynamicsFile))
    {
        luaerror(L, "cannot run config. file: %s", lua_tostring(L, -1));
        exit(EXIT_FAILURE);
    }


    int luaquit = 0;
    size_t iter = 0;
    do{
        iter++;
        p->iter_final++;

        /* Get settings from the lua script */
        lua_getglobal(L, "getConfig"); /* function to be called */
        // Push function arguments
        lua_pushnumber(L, iter);
        lua_pushnumber(L, p->newx);
        lua_pushnumber(L, p->n_beads);

        /* do the call (2 arguments, 0 result) */
        if (lua_pcall(L, 3, 0, 0) != LUA_OK)
        {
            luaerror(L, "error running function 'getConfig' in %s : %s",
                     p->luaDynamicsFile,
                     lua_tostring(L, -1));
            exit(EXIT_FAILURE);
        }

        /* retrieve result */
        fconf.kDom = lua_get_float(L, "kDom");
        fconf.kVol = lua_get_float(L, "kVol");
        fconf.kInt = lua_get_float(L, "kInt");
        fconf.kBackbone = lua_get_float(L, "kBackbone");
        fconf.kRad = lua_get_float(L, "kRad");
        fconf.kBeadWell = lua_get_float(L, "kBeadWell");
        fconf.kChrWell = lua_get_float(L, "kChrWell");
        fconf.top_plane = lua_get_float(L, "top_plane");
        fconf.bottom_plane = lua_get_float(L, "bottom_plane");
        p->compress = lua_get_float(L, "kCom");
        double Fb = lua_get_float(L, "fBrown");
        luaquit = lua_get_int(L, "quit");
        fconf.dInteraction = fconf.r0 * lua_get_float(L, "dInteraction");

        if(p->verbose > 10)
        {
            printf("kInteraction: %f\n", fconf.kInt);
            printf("kBackbone:    %f\n", fconf.kBackbone);
            printf("dInteraction: %f\n", fconf.dInteraction);
        }

        /* Start of Molecular Dynamics
         * 2.1 Gradient of functional */
        grad3(X,
              p->n_beads,
              p->R,
              p->I,
              p->active_pair,
              p->backbone,
              p->n_backbone,
              g,
              &fconf);

        /* Cap gradient for stability */
        for(size_t kk = 0; kk<3*p->n_beads; kk++)
        {
            if(fabs(g[kk]) > 1.0)
            {
                g[kk] = copysign(1.0, g[kk]);
            }
        }

        if(p->compress > 0)
        {
            comforce(p, g);
        }

        /* 2.2 Brownian force */
        if(Fb>0)
        {
            for(size_t kk = 0; kk < p->n_beads; kk++)
            {
                for(size_t idx = 0 ; idx<3; idx++)
                {
                    g[3*kk+idx] += 0.1*Fb*normal();
                }
            }
        }

        /*
         * Bead wells i.e. attraction of beads to specific coordinates
         *
         */
        if( p->n_bead_wells > 0 )
        {
            bead_wells_gradient(&fconf, p->n_beads, p->bead_wells, p->n_bead_wells, X, g);
        }

        /* 2.3 Dampening */
        /* Estimate velocities,
           note: cheating and doing it per component */

        if(damp > 0)
        {
            for(size_t pp = 0 ; pp < 3*p->n_beads ; pp++)
            {
                v[pp] = (X[pp] - Xm[pp]) / (2.0 * dt);
                g[pp] = g[pp] + damp*v[pp];
            }
        }

        /* 3. Update X */
        for(size_t pp = 0 ; pp < 3*p->n_beads ; pp++)
        {
            double xt = X[pp];
            X[pp] = 2.0*X[pp] - Xm[pp] - g[pp]*pow(dt,2);
            Xm[pp] = xt; // Update Xm to reflect the previous X-value
        }
        /* End of molecular dynamics */

        //
        // Enforce absolute bead placement if --absolute was provided
        //
        if(p->bead_apos != NULL)
        {
            for(i64 kk = 0; kk < p->n_bead_apos; kk++)
            {
                i64 bead_id = p->bead_apos[kk].bead_id;
                X[3*bead_id + 0] = p->bead_apos[kk].x;
                X[3*bead_id + 1] = p->bead_apos[kk].y;
                X[3*bead_id + 2] = p->bead_apos[kk].z;
            }
        }

        /*
         * Possibly output some info at the end of the step
         */
        if(iter % 1500 == 0 || iter == maxiter-1)
        { showStep = 1; } else { showStep = 0; }

        if(showStep == 1)
        {
            /* Calculate gradient 2-norm */
            gnorm = 0;
            for(size_t kk = 0; kk < 3*p->n_beads; kk++)
            {
                gnorm += pow(v[kk], 2);
            }
            gnorm = sqrt(gnorm);

            error = err3(X,
                         p->n_beads,
                         p->R,
                         p->I,
                         &fconf);
            mflock_logwrite(p, 2, "    Iter: %6zu, E: %e, ||G||: %e\n",
                            iter, error, gnorm);
            fflush(p->logf);
        }

    } while( (iter < maxiter) && (run == 1) && (luaquit == 0));

    // TODO this is already calculated at at end of the last iteration.
    // also DRY.

    /* At final step, report back */
    double errorFinal = err3(X,
                             p->n_beads,
                             p->R,
                             p->I,
                             &fconf);

    /* The gradient is a 3*p->n_beads-dimensional vector, we return the 2-norm
       as the grad_final */
    double gradFinal = 0;
    for(size_t kk = 0; kk < 3*p->n_beads; kk++)
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
    lua_close(L);

    return 0;
}

static void mflock_summary(mflock_t * p)
{
    const double * restrict X = p->beads;
    assert(p->n_beads>0);

    mflock_logwrite(p, 1, "\n");
    mflock_logwrite(p, 1, " >> Optimization summary:\n");
    mflock_logwrite(p, 1, "    Final iterations: %zu\n", p->iter_final);
    mflock_logwrite(p, 1, "    Total time: %zu s\n", p->time_final);
    mflock_logwrite(p, 1, "    Final gradient norm: %e\n", p->grad_final);
    mflock_logwrite(p, 1, "    Final total error: %e\n", p->err_final);
    size_t n_auto_used = 0;
    if(p->autocontacts){
        for(i64 kk = 0; kk < p->n_pairs; kk++) {
            if(p->active_pair[kk] == 1) {
                n_auto_used++;
            }
        }
    }
    mflock_logwrite(p, 1, "    Used %ld / %zu contact pairs\n", n_auto_used, p->n_pairs);

    // X: mean, max, min
    double mex = 0, mey = 0, mez = 0;
    double md = 10e99;
    double mix = md, miy = md, miz = md;    /* min, x, y, z */
    double max = -md, may = -md, maz = -md; /* max, x, y, z */
    double mer = 0, mir = 10e99, mar = 0;   /* mean, min, max of radius */

    for(size_t kk = 0 ; kk<p->n_beads; kk++)
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
        double pr = norm3d(X+3*kk);
        mar = dmax(mar, pr);
        mir = dmin(mir, pr);
        mer += pr;
    }
    mex /= p->n_beads;
    mey /= p->n_beads;
    mez /= p->n_beads;
    mer /= p->n_beads;

    mflock_logwrite(p, 2,  "\n");
    mflock_logwrite(p, 2,  ">> Structure summary:\n");
    mflock_logwrite(p, 2,  "              X       Y       Z       R\n");
    mflock_logwrite(p, 2,  "   Max:  % .3f, % .3f, % .3f, % .3f\n",
                    max, may, maz, mar);
    mflock_logwrite(p, 2,  "   Mean: % .3f, % .3f, % .3f, % .3f\n",
                    mex, mey, mez, mer);
    mflock_logwrite(p, 2,  "   Min:  % .3f, % .3f, % .3f, % .3f\n",
                    mix, miy, miz, mir);

    if(p->bead_wells != NULL)
    {
        i64 n_filled = 0;
        double r02 = pow(2.0*p->r0, 2.0);
        for(size_t kk = 0; kk < p->n_bead_wells; kk++)
        {
            int filled = 0;
            size_t accepts = p->bead_wells[kk].bead_idx;
            double * WX = (double*) &p->bead_wells[kk];
            double distance = eudist3p2(WX, X + 3*accepts);
#if 0
            printf("Well %zu (%f, %f, %f) accepts %zu: distance: %f\n", kk,
                   p->bead_wells[kk].P.x,
                   p->bead_wells[kk].P.y,
                   p->bead_wells[kk].P.z,
                   accepts, distance/(2.0*r02));
#endif
            if(distance < 2.0*r02)
            {
                filled = 1;
            }
            if(p->diploid == 1)
            {
                accepts += p->n_beads/2;
                distance = eudist3p2(WX, X + 3*accepts);
#if 0
                printf("     %zu (%f, %f, %f) accepts %zu: distance: %f\n", kk,
                       p->bead_wells[kk].P.x,
                       p->bead_wells[kk].P.y,
                       p->bead_wells[kk].P.z,
                       accepts, distance/(2.0*r02));
#endif
                if(distance < 2.0*r02)
                {
                    filled = 1;
                }
            }

            if(filled > 0)
            {
                n_filled++;
            }
        }
        mflock_logwrite(p, 1, "%ld / %ld wells are filled with the accepted bead\n",
                        n_filled, p->n_bead_wells);
    }

    if(run == 0)
    {
        mflock_logwrite(p, 0, "abnormal exit (Ctrl+c was pressed?)\n");
    }



    char * timestr = cf_timestr();
    mflock_logwrite(p, 2, "Finished at: %s\n", timestr);
    free(timestr);
    return;
}

/* Check that the list of contacts is sorted and that there are no duplicates */
static void
check_bead_contacts(const u32 * pairs, i64 npairs,
                    const int verbose, const i64 nbead)
{
    if(verbose > 1)
    {
        printf("Validating %ld contact constraints\n", npairs);
    }
    for(i64 kk = 0; kk+1 < npairs; kk++)
    {
        if(pairs[2*kk] > pairs[2*(kk+1)])
        {
            printf("Contact pairs not sorted\n");
            exit(EXIT_FAILURE);
        }
        if(pairs[2*kk] == pairs[2*(kk+1)])
        {
            if(pairs[2*kk + 1] == pairs[2*(kk+1)+1])
            {
                printf("Contact pairs contains a duplicate\n");
                exit(EXIT_FAILURE);
            }
            if(pairs[2*kk + 1] > pairs[2*(kk+1) + 1])
            {
                printf("Contact pairs not sorted\n");
                exit(EXIT_FAILURE);
            }
        }
    }
    for(i64 kk = 0; kk < 2*npairs; kk++)
    {
        if(pairs[kk] >= nbead)
        {
            printf("Error: A contact refers to bead %u but there are only %ld\n", pairs[kk], nbead);
            exit(EXIT_FAILURE);
        }
    }
    if(verbose > 1)
    {
        printf("Contact pairs seems ok\n");
    }
}

static void
mflock_read_contact_pairs(mflock_t * p)
{
    if(p->contact_pairs_file == NULL)
    {
        mflock_logwrite(p, 1, "Pairwise interactions -- no file specified\n");
        return;
    }
    mflock_logwrite(p, 1, "Reading pairwise interactions from %s\n",
                    p->contact_pairs_file);
    uint64_t nCP = 0;
    if(npy_extension(p->contact_pairs_file))
    {
        int npair = 0;
        p->I = load_bead_contacts_from_npy(p->contact_pairs_file,
                                           &npair);
        if(p->I == NULL)
        {
            fprintf(stderr, "Failed to read contact pairs from %s\n",
                    p->contact_pairs_file);
            exit(EXIT_FAILURE);
        }
        p->n_pairs = npair;
        nCP = npair;
    } else {

        p->I = contact_pairs_read(p->contact_pairs_file, &nCP);
        if(p->I == NULL)
        {
            fprintf(stderr, "%s/%d Failed to read contact pairs from %s\n",
                    __FILE__, __LINE__, p->contact_pairs_file);
            exit(EXIT_FAILURE);
        }
        p->n_pairs = nCP;
    }

    check_bead_contacts(p->I, p->n_pairs, p->verbose, p->n_beads);

    mflock_logwrite(p, 1, "Read %lu contacts pairs\n", nCP);

    return;
}

static void mflock_init_coordinates(mflock_t * p)
{
    if(p->verbose > 2)
    {
        printf("init coordinates\n");
        printf("   expecting %zu points\n", p->n_beads);
    }
    p->newx = 0;
    if(p->xfname != NULL)
    {
        p->beads = malloc(3*p->n_beads*sizeof(double));
        assert(p->beads != NULL);
        if(mflock_load_coordinates(p) != 0)
        {
            mflock_logwrite(p, 2, "Could not read coordinates from %s\n");
            exit(EXIT_FAILURE);
        }
        mflock_logwrite(p, 2, "Using coordinates from %s\n", p->xfname);
    }

    if(p->beads == NULL)
    {
        p->newx = 1;
        mflock_logwrite(p, 1, "Using random initialization for X\n");
        srand(p->rseed);
        p->beads = calloc(3*p->n_beads, sizeof(double));
        assert(p->beads != NULL);
        elli * E = p->E;
        if(E == NULL){
            E = elli_new(1, 1, 1);
        }

        for(size_t kk = 0; kk< p->n_beads; kk++)
        {
            int accepted = 0;
            while(accepted == 0)
            {
                for(int idx =0; idx<3; idx++)
                {
                    p->beads[3*kk+idx] = 2.0*(rand()/(double) RAND_MAX-.5);
                }
                if(elli_getScale(E, p->beads+3*kk) < 0.95)
                { accepted = 1;}
            }
        }

        if(p->E == NULL)
        {
            free(E);
        }
        mflock_logwrite(p, 1, "X[0] = %f\n", p->beads[0]);
    }
    if(p->verbose > 2)
    {
        printf("   coordinates loaded\n");
    }

    if(p->cmm_cmap != NULL)
    {
        mflock_logwrite(p, 2, "Loading color map from %s\n", p->cmm_cmap);
        p->cmap = load_cmap(p->cmm_cmap);
    }

    return;
}

static void mflock_validate_labels(const mflock_t * p)
{
    // Check how many polymers there are
    int npoly = 1;
    for(size_t kk = 1; kk< p->n_beads; kk++)
    {
        if(p->L[kk-1] != p->L[kk])
        {
            npoly++;
        }
    }


    mflock_logwrite(p, 1,
                    "%d polymers found in the labels\n", npoly);


    if( npoly > 55 )
    {
        if(p->verbose > 0)
        {
            printf("Found %d polymers / %d places where L[kk] != L[kk+1]. This probably "
                   "indicates that the label array is not constructed as expected by the "
                   "program.",
                   npoly, npoly-1);
        }
    }
    return;
}


static int mflock_load_bead_labels(mflock_t * p)
{
    if(p->lfname == NULL)
    {
        mflock_logwrite(p, 1, "No L-file specified\n");
        return EXIT_FAILURE;
    }

    mflock_logwrite(p, 1, "Reading L-labels from %s\n", p->lfname);

    if(npy_extension(p->lfname))
    {
        int nbead = 0;
        p->L = load_bead_labels_from_npy(p->lfname, &nbead);
        if(p->L == NULL)
        {
            fprintf(stderr, "Failed to read labels from %s\n\n",
                    p->lfname);
            exit(EXIT_FAILURE);
        }
        p->n_beads = nbead;
    } else {
        size_t fsize = cf_file_size(p->lfname);

        mflock_logwrite(p, 1, "As uint8_t, %s %zu numbers (%zu bytes)\n", p->lfname,
                        fsize/sizeof(uint8_t), fsize);

        // Try to read as binary
        FILE * f = fopen(p->lfname, "rb");
        if(f == NULL)
        {
            fprintf(stderr, "Unable to read %s\n", p->lfname);
            exit(EXIT_FAILURE);
        }

        p->n_beads = fsize/sizeof(uint8_t);
        if(p->n_beads == 0)
        {
            fprintf(stderr, "No beads desribed in %s\n", p->lfname);
            exit(EXIT_FAILURE);
        }

        p->L = malloc(p->n_beads*sizeof(double));
        assert(p->L != NULL);
        size_t nread = fread(p->L, sizeof(uint8_t), p->n_beads, f);
        if(nread != p->n_beads)
        {
            fprintf(stderr, "Unable to read from %s\n", p->lfname);
            exit(EXIT_FAILURE);
        }
        fclose(f);
    }

    mflock_validate_labels(p);

    if(p->diploid)
    {
        mflock_logwrite(p, 2, "Duplicating the beads for diploid structures\n");

        p->L = realloc(p->L, 2*p->n_beads*sizeof(uint8_t));
        assert(p->L != NULL);

        for(size_t kk = 0; kk < p->n_beads; kk++)
        {
            p->L[kk+p->n_beads] = p->L[kk] + 32;
        }
        p->n_beads = 2*p->n_beads;
    }

    mflock_logwrite(p, 2, "L = [%u, %u, ..., %u]\n", p->L[0], p->L[1], p->L[p->n_beads-1]);

    return 0;
}


static int mflock_load_radial_constraints(mflock_t * p)
{

    if(p->rfname == NULL)
    {
        mflock_logwrite(p, 1, "No radial preferences to read\n");
        return EXIT_FAILURE;
    }

    mflock_logwrite(p, 1, "Reading R-values from %s\n", p->rfname);
    size_t nbytes = 0;
    p->R = (double *) wio_read(p->rfname, &nbytes);
    printf("Read %zu doubles from %s\n", nbytes/sizeof(double), p->rfname);

    if(2*nbytes/sizeof(double) == p->n_beads)
    {
        printf("Found half as many values as beads, assuming diploid and duplicating data\n");
        if(p->n_beads%2 == 1)
        {
            printf("ERROR: N%%2 == 1\n");
            exit(-1);
        }

        p->R = realloc(p->R, p->n_beads*sizeof(double));
        if(p->R == NULL)
        {
            printf("ERROR: out of memory\n");
            exit(-1);
        }
        for(size_t pp = 0; pp < p->n_beads/2; pp++)
        {
            p->R[pp + p->n_beads/2] = p->R[pp];
        }
        //memcpy(p->R+p->n_beads/2*sizeof(double), p->R, p->n_beads/2*sizeof(double));
    }


    if(!( (nbytes/sizeof(double) == p->n_beads) || (nbytes/sizeof(double) == p->n_beads/2) ))
    {
        printf("Error: Can't make sense of %s, expected %zu or %zu bytes but got %zu\n", p->rfname, 4*p->n_beads, 2*p->n_beads, nbytes);
        exit(-1);
    }

    size_t nInf = 0;
    for(size_t kk = 0; kk<p->n_beads; kk++)
    {
        if(!isfinite(p->R[kk]))
        {
            nInf++;
        }
    }

    printf("%s contains %zu non-finite values which will be ignored.\n", p->rfname, nInf);

    return EXIT_SUCCESS;
}



static void mflock_show(mflock_t * p, FILE * f)
{
    double volocc = p->n_beads*4.0/3.0*M_PI*pow(p->r0,3) / (4.0/3.0*M_PI);

    if(p->E != NULL)
    {
        volocc = p->n_beads*4.0/3.0*M_PI*pow(p->r0,3) / elli_vol(p->E);
    }

    fprintf(f, "\n");
    fprintf(f, " >> Parameters:\n");
    fprintf(f, "    Problem size: %zu points (%zu variables)\n", p->n_beads, 3*p->n_beads);
    fprintf(f, "    Dynamics program: %s\n", p->luaDynamicsFile);
    fprintf(f, "    Bead radius %f (vol. occ. %f)\n", p->r0, volocc);
    fprintf(f, "    random seed: %zu\n", p->rseed);
    fprintf(f, "    write compresseed cmm: %d\n", p->cmmz);

    fprintf(f, "    geometry: ");
    switch(p->geometry)
    {
    case MFLOCK_BOX:
        fprintf(f, "Box [-1, 1]^3\n");
        break;
    case MFLOCK_SPHERE:
        fprintf(f, "Sphere; radius = 1, center = (0, 0, 0))\n");
        break;
    case MFLOCK_ELLIPSOID:
        assert(p->E != NULL);
        fprintf(f, "Ellipsoid; axes = %f, %f, %f, center =  (0, 0, 0)\n",
                p->E->a, p->E->b, p->E->c);
        break;
    }
    if(p->contact_pairs_file == NULL)
    {
        fprintf(f, "    Contact Pairs file not specified\n");
    } else {
        fprintf(f, "   Contacts Pairs: %s\n", p->contact_pairs_file);
    }

    if(p->xfname == NULL)
    {
        fprintf(f, "    Coordinates not loaded, will use random initialization\n");
    }
    else
    {
        fprintf(f, "    Coordinates from: %s\n", p->xfname);
    }

    if(p->lfname == NULL)
    {
        fprintf(f, "    Labels not provided\n");
    }
    else
    {
        fprintf(f, "    Labels from: %s\n", p->lfname);
    }

    if(p->rfname == NULL)
    {
        fprintf(f, "    Radial preferences not set\n");
    }
    else
    {
        fprintf(f, "    Radial preferences: %s\n", p->rfname);
    }

    if(p->fname_bead_wells == NULL)
    {
        fprintf(f, "    Bead wells: -not specified-\n");
    }
    else
    {
        fprintf(f, "    Bead wells: %s\n", p->fname_bead_wells);
    }

    if(p->ofoldername == NULL)
    {
        fprintf(f, "    output folder not specified\n");
    }
    else
    {
        fprintf(f, "    output folder: %s\n", p->ofoldername);
    }

    fprintf(f, "    diploid=%s\n", cf_YES_NO(p->diploid));

    fprintf(f, "    live view=%s\n", cf_YES_NO(p->liveView));

    if(p->liveView) {
        fprintf(f, "    auto close=%s\n", cf_YES_NO(p->live_auto_close));
    }

    fprintf(f, "    verbose level: %d\n", p->verbose);
    fprintf(f, "\n");
    return;
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
static void dump_lua_dynamics(const char * luafile)
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

            double kDom = lua_get_float(L, "kDom");
            double kVol = lua_get_float(L, "kVol");
            double dInt = lua_get_float(L, "dInteraction");
            double kInt = lua_get_float(L, "kInt");
            double kRad = lua_get_float(L, "kRad");
            double compress = lua_get_float(L, "kCom");
            double Fb = lua_get_float(L, "fBrown");
            luaquit = lua_get_int(L, "quit");

            printf("%d, %zu, %f, %f, %f, %f, %f, %f, %f\n",
                   newx, iter, kDom, kVol, kInt, dInt, kRad, compress, Fb);
        } while(luaquit == 0);

        lua_close(L);
    }
    return;
}
static void test_read_write_csv(void)
{
    int nbead = 71;
    double * X = calloc(3*nbead, sizeof(double));
    double * Y = calloc(3*nbead, sizeof(double));
    for(i64 kk = 0; kk < 3*nbead; kk++)
    {
        X[kk] = 2.0*((double) rand() / (double) RAND_MAX) - 1.0;
    }

    char * tmpfile = calloc(1024, 1);
    snprintf(tmpfile, 1024, "tmp_XXXXXX");
    int fd = mkstemp(tmpfile);
    if(fd == -1)
    {
        fprintf(stderr, "Failed to create a temporary file\n");
        exit(EXIT_FAILURE);
    }
    close(fd);

    if(write_bead_coordinates_to_csv(tmpfile, X, nbead, MFLOCK_SPHERE, NULL))
    {
        exit(EXIT_FAILURE);
    }

    if(load_bead_coordinates_from_csv(tmpfile, NULL, Y, nbead))
    {
        exit(EXIT_FAILURE);
    }

    int err = 0;
    for(int kk = 0; kk < nbead; kk++)
    {
        for(int ii = 0; ii < 3; ii++)
        {
            if( fabs(X[3*kk + ii] - Y[3*kk + ii]) > 1e-4)
            { err = 1; }
        }
    }
    remove(tmpfile);

    if(err)
    {
        fprintf(stderr, "Failed\n"
                "write_bead_coordinates_to_csv(tmpfile, X, 1, NULL);\n"
                "load_bead_coordinates_from_csv(tmpfile, Y, 1);\n");
        exit(EXIT_FAILURE);
    }
    free(X);
    free(Y);
    free(tmpfile);
    return;
}

static void mflock_ut(void)
{
    printf("Testing ... \n");
    int nfail = 0;
    if(npy_extension("a.npy") == 0)
    {
        nfail++;
    }
    if(npy_extension("abc.NPY") == 0)
    {
        nfail++;
    }
    if(npy_extension("abc.PY") == 1)
    {
        nfail++;
    }

    test_read_write_csv();

    printf("All tests passed\n");
    return;
}

static void
mflock_usage()
{
    printf("mflock %s ", cf_version);
    printf("%s\n", src_mflock_help_txt);
    return;
}


static mflock_cli_status
mflock_parse_cli(mflock_t * p, int argc, char ** argv)
{
    /* Specifications of ellipsoid */
    double ea = -1;
    double eb = -1;
    double ec = -1;

    struct option longopts[] = {
        { "version",       no_argument,       NULL,   'i' },
        { "help",          no_argument,       NULL,   'h' },
        { "test",          no_argument,       NULL,   'T' },
        /* Data */
        { "contact-pairs", required_argument, NULL,   'p' },
        { "autocontacts",  no_argument,       NULL,   '2' },
        { "xFile",         required_argument, NULL,   'x' },
        { "coordinates",   required_argument, NULL,   'x' },
        { "rFile",         required_argument, NULL,   'r' },
        { "radii",         required_argument, NULL,   'r' },
        { "labels",        required_argument, NULL,   'L' },
        { "lFile",         required_argument, NULL,   'L' },
        { "outFolder",     required_argument, NULL,   'o' },
        { "bead-wells",    required_argument, NULL,   'W' },
        { "absolute",      required_argument, NULL,   'P' },
        /* Settings */
        { "backbone",      no_argument,       NULL,   'b' },
        { "diploid",       no_argument,       NULL,   'D' },
        { "maxiter",       required_argument, NULL,   'n' },
        { "maxtime",       required_argument, NULL,   't' },
        { "seed",          required_argument, NULL,   's' },
        { "verbose",       required_argument, NULL,   'v' },
        { "live",          no_argument,       NULL,   'a' },
        { "liveclose",     no_argument,       NULL,   'e' },
        { "cmm",           no_argument,       NULL,   'c' },
        { "cmmz",          no_argument,       NULL,   'z' },
        { "cmap",          required_argument, NULL,   '1' },
        { "defaults",      no_argument,       NULL,   'd' },
        /* Geometry */
        { "box",           no_argument,       NULL,   'X' },
        { "radius",        required_argument, NULL,   'R' },
        { "vq",            required_argument, NULL,   'Q' },
        { "ea",            required_argument, NULL,   'A' },
        { "eb",            required_argument, NULL,   'B' },
        { "ec",            required_argument, NULL,   'C' },

        // Lua program for dynamics
        { "dconf",         required_argument, NULL,   'l' },
        { "config",        required_argument, NULL,   'l' },
        { "dconf-show",    required_argument, NULL,   'M' },
        { "config-show",   required_argument, NULL,   'M' },

        // debug
        { "use-csv",       no_argument,       NULL,    'u' },
        { NULL,            0,                 NULL,    0  }
    };

    int ch;
    while((ch = getopt_long(argc, argv,
                            "12:abA:B:cC:De:w:x:r:n:p:P:t:R:v:o:hMs:L:zcdQ:l:W:TuX",
                            longopts, NULL)) != -1)
    {
        switch(ch) {
        case '1':
            free(p->cmm_cmap);
            p->cmm_cmap = strdup(optarg);
            break;
        case '2':
            p->autocontacts = 1;
            break;
        case 'a':
            p->liveView = 1;
            break;
        case 'b':
            p->create_backbone = 1;
            break;
        case 'A':
            ea = atof(optarg);
            p->geometry = MFLOCK_ELLIPSOID;
            break;
        case 'B':
            eb = atof(optarg);
            p->geometry = MFLOCK_ELLIPSOID;
            break;
        case 'c':
            p->write_cmm = 1;
            break;
        case 'C':
            ec = atof(optarg);
            p->geometry = MFLOCK_ELLIPSOID;
            break;
        case 'i':
            printf("mflock (chromflock version %s)\n", cf_version);
            printf("Build date: %s, %s\n", __DATE__, __TIME__);
            return MFLOCK_ARGS_QUIT;
        case 'd':
            printf("Defaults:\n");
            mflock_show(p, stdout);
            return MFLOCK_ARGS_QUIT;
        case 'D':
            p->diploid = 1;
            break;
        case 'e':
            p->liveView = 1;
            p->live_auto_close = 1;
            break;
        case 'W':
            free(p->fname_bead_wells);
            p->fname_bead_wells = strdup(optarg);
            assert(p->fname_bead_wells != NULL);
            break;
        case 'x':
            free(p->xfname);
            p->xfname = strdup(optarg);
            assert(p->xfname != NULL);
            break;
        case 'r':
            free(p->rfname);
            p->rfname = strdup(optarg);
            assert(p->rfname != NULL);
            break;
        case 'L':
            free(p->lfname);
            p->lfname = strdup(optarg);
            assert(p->lfname != NULL);
            break;
        case 'n':
            p->maxiter = atol(optarg);
            break;
        case 'p':
            free(p->contact_pairs_file);
            p->contact_pairs_file = strdup(optarg);
            break;
        case 'P':
            free(p->bead_apos_file);
            p->bead_apos_file = strdup(optarg);
            break;
        case 't':
            p->maxtime = atol(optarg);
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
        case 'u':
            p->use_csv = 1;
            break;
        case 'v':
            p->verbose = atoi(optarg);
            break;
        case 'o':
            free(p->ofoldername);
            p->ofoldername = strdup(optarg);
            assert(p->ofoldername != NULL);
            break;
        case 'l':
            free(p->luaDynamicsFile);
            p->luaDynamicsFile = strdup(optarg);
            assert(p->luaDynamicsFile != NULL);
            break;
        case 'M':
            dump_lua_dynamics(optarg);
            return MFLOCK_ARGS_QUIT;
        case 'h':
            return(1);
        case 'X':
            p->geometry = MFLOCK_BOX;
            break;
        case 'z':
            p->write_cmm = 1;
            p->cmmz = 1;
            break;
        case 'T':
            mflock_ut();
            return MFLOCK_ARGS_QUIT;
        default:
            return MFLOCK_ARGS_ERR;
        }
    }

    if(p->bead_apos_file != NULL)
    {
        FILE * fid = fopen(p->bead_apos_file, "rb");
        if(fid == NULL)
        {
            printf("Unable to open %s\n", p->bead_apos_file);
            return MFLOCK_ARGS_ERR;
        }
        fclose(fid);
    }

    if(p->luaDynamicsFile == NULL)
    {
        fprintf(stderr, "Incorrect command line: --dconf not specified\n");
        fprintf(stderr, "see %s --help for usage\n", argv[0]);
        exit(EXIT_FAILURE);
    }

    /* Make sure that the outfolder ends with a path_separator.
       Note: it should already be allocated to have room for it if missing */
    if(p->ofoldername)
    {
        if(p->ofoldername[strlen(p->ofoldername) - 1] != '/' )
        {
            p->ofoldername = realloc(p->ofoldername, strlen(p->ofoldername)+1);
            p->ofoldername[strlen(p->ofoldername)+1] = '\0';
            p->ofoldername[strlen(p->ofoldername)] = '/';
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
        printf("ERROR: Ellipsoidal geometry requires 1 = ea >= eb >= ec > 0"
               " \n\n");
        return MFLOCK_ARGS_ERR;
    }

    return MFLOCK_ARGS_OK;
}

static mflock_t * mflock_new(void)
{
    mflock_t * p = calloc(1, sizeof(mflock_t));
    assert(p != NULL);

    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);

    // TODO: overflows
    p->rseed = time(NULL)*getpid()*ts.tv_nsec;

    p->maxiter = 1000000; // iterations
    p->maxtime = 60*60*10; // seconds

    p->r0 = -1;
    p->volq = 0.2;

    p->verbose = 1;
    p->newx = 1;
    p->write_cmm = 0;
    p->geometry = MFLOCK_SPHERE;

    /* Create a suggestion for the output folder */
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

static void mflock_set_and_create_output_folder(mflock_t * mf)
{
    struct stat st;
    memset(&st, 0, sizeof(struct stat));

    if(stat(mf->ofoldername, &st) == -1)
    {
        if(mf->verbose >= 1){
            printf("Creating output folder: %s\n", mf->ofoldername);
        }
        int folderok = mkdir(mf->ofoldername, 0770);
        if(folderok != 0)
        {
            printf("Could not create output folder\n");
            exit(-1);
        }

    } else {
        if(mf->verbose >= 2){
            printf("Output folder did already exist.\n");
        }
    }

    char * fullofolder = realpath(mf->ofoldername, NULL);
    if(fullofolder == NULL)
    {
        printf("Error: The output folder '%s' can't be understood by realpath\n", mf->ofoldername);
        exit(-1);
    }

    // Create output folder
    if(mf->verbose >= 1) {
        printf("Output folder: %s\n", fullofolder);
    }
    free(fullofolder);
    return;
}

static void mflock_load_bead_wells(mflock_t * mf)
{
    // side effects:
    // mf->n_bead_wells
    // mf->bead_wells

    if(mf->fname_bead_wells == NULL)
    {
        if(mf->verbose > 1)
        {
            printf("no bead wells to load\n");
        }
        return;
    }

    if(mf->verbose > 1)
    {
        printf("loading bead wells from %s\n", mf->fname_bead_wells);
    }

    if(npy_extension(mf->fname_bead_wells))
    {
        int n_bead_wells = 0;
        bpos * pos = load_bead_apos_from_npy(mf->fname_bead_wells, &n_bead_wells);

        if(pos == NULL)
        {
            fprintf(stderr, "Failed to load bead wells from %s\n",
                    mf->fname_bead_wells);
            exit(EXIT_FAILURE);
        }

        mf->bead_wells = calloc(n_bead_wells, sizeof(wpos));
        assert(mf->bead_wells != NULL);
        mf->n_bead_wells = n_bead_wells;
        wpos * wells = (wpos *) mf->bead_wells;
        for(i64 kk = 0; kk < n_bead_wells; kk++)
        {
            wells[kk].bead_idx = (size_t) pos[kk].bead_id;
            wells[kk].P.x = (double) pos[kk].x;
            wells[kk].P.y = (double) pos[kk].y;
            wells[kk].P.z = (double) pos[kk].z;
        }
        free(pos);
        return;
    }


    size_t fsize = cf_file_size(mf->fname_bead_wells);
    if(fsize % sizeof(double) != 0)
    {
        fprintf(stderr, "%s seems corrupt, The file size (%zu) can not be divided by %zu",
                mf->fname_bead_wells, fsize, sizeof(double));
        exit(EXIT_FAILURE);
    }

    size_t n_elements = fsize / sizeof(double);

    mf->bead_wells = calloc(n_elements, sizeof(double));
    FILE * fid = fopen(mf->fname_bead_wells, "r");
    if(fid == NULL)
    {
        fprintf(stderr, "Error while reading %s\n", mf->fname_bead_wells);
        exit(EXIT_FAILURE);
    }
    size_t n_read = fread(mf->bead_wells, sizeof(double), n_elements, fid);
    if(n_read != n_elements)
    {
        fprintf(stderr, "Error while reading %s\n", mf->fname_bead_wells);
        exit(EXIT_FAILURE);
    }
    fclose(fid);

    size_t n_constraints = n_elements / 4;
    for(size_t kk = 0; kk < n_constraints ; kk++)
    {
        size_t bead1 = mf->bead_wells[kk].bead_idx;
        if(bead1 > mf->n_beads)
        {
            fprintf(stderr, "Error: Got a bead well for bead %zu, but there are only %zu beads\n",
                    bead1, mf->n_beads);
            exit(EXIT_FAILURE);
        }
    }
    mf->n_bead_wells = n_constraints;
    if(mf->verbose > 1)
    {
        printf("Successfully read %u beads wells from %s\n",
               mf->n_bead_wells,
               mf->fname_bead_wells);
    }
    return;
}

static void
mflock_update_autocontacts(mflock_t * mf)
{
    double dist2_th = pow(mf->r0*3.0, 2.0);
    int n_activated = 0;
    for(size_t pp = 0; pp < mf->n_pairs; pp++)
    {
        u32 * pair = mf->I + 2*pp;
        double dist2 = eudist3p2(mf->beads + 3*pair[0], mf->beads + 3*pair[1]);

        if(dist2 < dist2_th)
        {
            mf->active_pair[pp] = 1;
            n_activated++;
        }
    }
    mflock_logwrite(mf, 2, "Activated %d / %zu contact pairs\n", n_activated, mf->n_pairs);
    return;
}

static void
mflock_init_autopairs(mflock_t * mf)
{
    if(mf->autocontacts == 0) {
        mflock_logwrite(mf, 2, "No autocontacts\n");
        return; // unwanted
    }

    mflock_logwrite(mf, 2, "Initializing autocontacts\n");
    mf->active_pair = calloc(mf->n_pairs, sizeof(u8));

    if(mf->newx == 0) {
        mflock_update_autocontacts(mf);
    }
    return;
}

// If mf->autoconfig, this will write down an array containing
// [[enforced, distance], [enforced, distance] ... ] where enforced
// is either 1 if the the contact was in use at the end of the run, or
// 0. Distance is the distance between the two corresponding beads
static void
mflock_write_autopairs(mflock_t * mf)
{
    if(mf->autocontacts == 0) {
        return;
    }
    assert(mf->active_pair != NULL);

    f32 * contact_details = malloc(mf->n_pairs*2*sizeof(u32));
    for(i32 kk = 0; kk < mf->n_pairs; kk++)
    {
        contact_details[2*kk] = mf->active_pair[kk];
        size_t a = mf->I[2*kk];
        size_t b = mf->I[2*kk+1];
        contact_details[2*kk+1] = sqrt(eudist3p2(mf->beads + 3*a , mf->beads + 3*b));
    }

    // TODO: npio_write
    char * fname = malloc(1024);
    sprintf(fname, "%s%s", mf->ofoldername, "contact_details.npy");
    int shape[2] = {mf->n_pairs, 2};
    i64 status = npio_write(fname, 2,
                            shape,
                            (const void *) contact_details,
                            NPIO_F32, NPIO_F32);
    assert(status != -1);
    free(fname);
    free(contact_details);
}


static void mflock_init_backbone(mflock_t * mf)
{
    if(mf->create_backbone == 0) {
        return;
    }

    mf->backbone = malloc(2*mf->n_beads*sizeof(u32));

    int n_backbone = 0;
    for(size_t kk = 0; kk + 1 < mf->n_beads; kk++)
    {
        if(mf->L[kk] == mf->L[kk+1])
        {
            mf->backbone[2*n_backbone + 0] = kk;
            mf->backbone[2*n_backbone + 1] = kk + 1;
            n_backbone++;
        }
    }
    mf->n_backbone = n_backbone;
    mflock_logwrite(mf, 2, "Added %d backbone contacts\n", mf->n_backbone);
    return;
 }

static void mflock_init(mflock_t * mf, int argc, char ** argv)
{
    mflock_set_and_create_output_folder(mf);

    /* Set names of output files */
    mf->xoutfname = malloc(1024*sizeof(char));
    assert(mf->xoutfname != NULL);
    if(mf->use_csv)
    {
        sprintf(mf->xoutfname, "%s%s", mf->ofoldername, "coords.csv");
    } else {
        sprintf(mf->xoutfname, "%s%s", mf->ofoldername, "coords.npy");
    }

    mf->logfname = malloc(1024*sizeof(char));
    assert(mf->logfname != NULL);
    sprintf(mf->logfname, "%s%s", mf->ofoldername, "log.txt");

    /* Open log */
    mf->logf = fopen(mf->logfname, "a");
    assert(mf->logf != NULL);

    if(mf->logf== NULL)
    {
        fprintf(stderr, "mflock: Failed to open log file for writing (%s)\n",
                mf->logfname);
        exit(EXIT_FAILURE);
    }

    char * time_str = cf_timestr();
    fprintf(mf->logf, "\nmflock started: %s\n", time_str);
    free(time_str);

    fprintf(mf->logf, "CMD: ");
    for(int kk = 0; kk<argc; kk++)
    {
        fprintf(mf->logf, "'%s' ", argv[kk]);
    }
    fprintf(mf->logf, "\n");
    fflush(mf->logf);

    /* We load the labels first to determine how many beads there
     * are (or 2X if --diploid is set) */
    mflock_load_bead_labels(mf);
    mflock_logwrite(mf, 1, "n_beads = %zu\n", mf->n_beads);


    /* Once we know how many beads we can set their radius based on the
     * volume quotient (or do nothing if it was given at the command line) */
    mflock_set_bead_size(mf);

    /* Read coordinates from previous iteration, or generate
       random coordinates if first iteration.
    */
    mflock_init_coordinates(mf);

    mflock_read_contact_pairs(mf);

    // Activate pairs based on distance in current
    // structure
    // TODO: Split to init and update
    mflock_init_autopairs(mf);

    // Add contacts between adjacent beads (--backbone)
    mflock_init_backbone(mf);

    mflock_load_radial_constraints(mf);

    mflock_load_bead_wells(mf);

    if(mf->bead_apos_file != NULL)
    {
        // bpos*
        mf->bead_apos = mflock_load_bead_apos(mf->bead_apos_file,
                                              &mf->n_bead_apos);

        if(mf->bead_apos == NULL)
        {
            fprintf(stderr,
                    "Unable to load absolut bead positions from %s\n",
                    mf->bead_apos_file);
            exit(EXIT_FAILURE);
        }
        // Check that the apos refers to existing beads
        for(i64 kk = 0; kk < mf->n_bead_apos; kk++)
        {
            if(mf->bead_apos[kk].bead_id >= mf->n_beads)
            {
                fprintf(stderr,
                        "Error: absolute position %ld refers to bead %d, "
                        "but there are only %ld beads\n",
                        kk+1, // one indexed
                        mf->bead_apos[kk].bead_id+1, // one-indexed
                        mf->n_beads);
                exit(EXIT_FAILURE);
            }
        }
    }


    if(mf->verbose > 3)
    {
        // show the loaded beads
        for(size_t kk = 0; kk< mf->n_beads; kk++)
        {
            printf("%5zu: %f, %f, %f\n",
                   kk, mf->beads[3*kk], mf->beads[3*kk+1], mf->beads[3*kk+2]);
        }
        // show loaded contacts
        for(size_t kk = 0; kk < mf->n_pairs; kk++)
        {
            printf("P%5zu: %u %u\n", kk,
                   mf->I[2*kk], mf->I[2*kk+1]);
            assert(mf->I[2*kk] < mf->n_beads);
            assert(mf->I[2*kk+1] < mf->n_beads);
        }
    }

    return;
}


void mflock_free(mflock_t * p)
{
    free(p->R);
    free(p->I);
    free(p->backbone);
    free(p->L);
    free(p->contact_pairs_file);
    free(p->lfname);
    free(p->rfname);
    free(p->xfname);
    free(p->xoutfname);
    free(p->ofoldername);
    free(p->luaDynamicsFile);
    free(p->logfname);
    free(p->E);
    free(p->beads);
    free(p->fname_bead_wells);
    free(p->bead_wells);
    free(p->bead_apos);
    free(p->bead_apos_file);
    free(p->cmm_cmap);
    free(p->cmap);
    free(p->active_pair);
    free(p);
    return;
}


static void
mflock_set_bead_size(mflock_t * p)
{
    if(p->r0 > 0) {
        return;
    }

    double vq = p->volq;

    mflock_logwrite(p, 2,
                    "Bead radius not given, setting total bead volume to %.2f%%\n",
                    100*vq);

    double Vd = 4.0/3.0*M_PI;
    if(p->E != NULL)
    {
        Vd = elli_vol(p->E);
    }
    p->r0 = cbrt( 3.0*vq*Vd / (4.0*p->n_beads*M_PI) );


    double Vbeads = p->n_beads*pow(p->r0,3)*M_PI*4.0/3.0;
    mflock_logwrite(p, 2,
                    "Vd: %f, Vbeads: %f, Vbeads/Vd = %f\n",
                    Vd, Vbeads, Vbeads/Vd);
    assert(fabs(Vbeads/Vd - vq)<1e-5);
    return;
}

// Write both to std out and to the log file (p->logfile)
// Always writes to the log file
// Write to stdout if p->verbose >= level
// When level == 0, a "WARNING: " will be pre-pended to the message
static void
mflock_logwrite(const mflock_t * p,
                int level,
                const char *fmt, ...)
{
    va_list args, args2;
    va_start(args, fmt);
    va_copy(args2, args);

    // Write to terminal
    if(p->verbose >= level) {
        fprintf(stdout, "    ");

        if(level == 0) {
            fprintf(stdout, "WARNING: ");
        }
        vfprintf(stdout, fmt, args);
    }

    // Write to log file
    if(p->logf != NULL) {
        if(level == 0) {
            fprintf(p->logf, "WARNING: ");
        }
        vfprintf(p->logf, fmt, args2);
    }
    va_end(args2);
    va_end(args);
    return;
}

// Load an expected number of beads into a pre-allocated memory location
static int mflock_load_coordinates(mflock_t * p)
{
    if(npy_extension(p->xfname))
    {
        return load_bead_coordinates_from_npy(p->xfname,
                                              NULL, p->beads,
                                              p->n_beads);
    } else {
        return load_bead_coordinates_from_csv(p->xfname,
                                              NULL, p->beads,
                                              p->n_beads);
    }
}

static int mflock_save_coordinates(mflock_t * p)
{
    if(run == 1)
    {
        if(p->verbose > 1)
        {
            mflock_logwrite(p, 1, "Writing final structure to %s\n", p->xoutfname);
        }
    }
    else {
        p->xoutfname = realloc(p->xoutfname, sizeof(char)*(strlen(p->xoutfname)+5));
        strcat(p->xoutfname, ".0");
        //    sprintf(p->xoutfname, "%s.0", p->xoutfname);

        mflock_logwrite(p, 1, "Writing non-finished structure to: %s\n", p->xoutfname);
    }

    if(p->verbose > 1)
    {
        mflock_logwrite(p, 1, "Columns: x, y, z, r\n");
    }


    if(p->use_csv)
    {
        return write_bead_coordinates_to_csv(p->xoutfname,
                                             p->beads,
                                             p->n_beads,
                                             p->geometry,
                                             p->E);
    } else {
        return write_bead_coordinates_to_npy(p->xoutfname,
                                             p->beads,
                                             p->n_beads,
                                             p->geometry,
                                             p->E);
    }

}


/* Wrapper if starting the main loop from a thread */
static void * solve_t(void * args)
{
    mflock_t * p = (mflock_t *) args;
    mflock_dynamics(p);
    if(p->live_auto_close) {
        p->quit_live_view = 1;
    }
    return NULL;
}

/* Write chimera file (for simple visualization) */
static void mflock_write_cmm(const mflock_t * p)
{
    if(p->write_cmm == 0)
    {
        return;
    }
    const double * restrict X = p->beads;

    if(p->cmmz == 1)
    {
        char * cmmfile = malloc(1024*sizeof(char));
        assert(cmmfile != NULL);
        sprintf(cmmfile, "%s/cmmdump.cmm.gz", p->ofoldername);

        cmmwritez(cmmfile, X, p->n_beads, p->r0, p->backbone, p->n_backbone, p->L, p->cmap);
        free(cmmfile);
    } else {
        char * cmmfile = malloc(1024*sizeof(char));
        assert(cmmfile != NULL);
        sprintf(cmmfile, "%s/cmmdump.cmm", p->ofoldername);
        cmmwrite(cmmfile, X, p->n_beads, p->r0, p->backbone, p->n_backbone, p->L, p->cmap);
        free(cmmfile);
    }

    return;
}

static void mflock_close_log(mflock_t * p)
{


    char * time_str = cf_timestr();
    fprintf(p->logf, "\nmflock finished: %s\n", time_str);
    free(time_str);

    /* Close log file */
    fclose(p->logf);
    return;
}

static void mflock_run(mflock_t * p)
{

    time_t starttime, nowtime;
    time(&starttime);

    /* Set up capturing of Ctrl+C for graceful exit */
    struct sigaction act;
    memset (&act, '\0', sizeof(act));
    act.sa_handler = stoprun;
    sigaction(SIGINT, &act, NULL);

    /* Start molecular dynamics */
    mflock_logwrite(p, 1, " >> Solving ... \n");


#ifdef SDL
    if(p->liveView == 1)
    {


        /* SDL likes to be in the main thread so we run mflock
           dynamics in a secondary.  */

        pthread_t th;

        pthread_create(&th, // thread
                       NULL, // pthread_attrib_t
                       solve_t, // function
                       p); // arg

        elli * E;
        if(p->E == NULL)
        {
            E = elli_new(1,1,1);
        } else {
            E = p->E;
        }

        liveview(p->beads, p->L, p->n_beads,
                 &p->quit_live_view, p->r0, E, p->cmap);
        pthread_join(th, NULL);
    } else {
        mflock_dynamics(p);
    }
#else
    if(p->liveView == 1) {
        printf("WARNING: Can't open the live view (--live) since the program was not "
               "linked to SDL2\n\n");
    }
    mflock_dynamics(p);
#endif


    /* resolution: 1 s, consider clockdiff  */
    time(&nowtime);
    p->time_final = difftime(nowtime, starttime);
    return;
}

int mflock(int argc, char ** argv)
{

#ifndef NDEBUG
    printf("WARNING: mflock was compiled without defining NDEBUG and will "
           "be very slow\n");
#endif

    mflock_t * mf = mflock_new();

    switch(mflock_parse_cli(mf, argc, argv))
    {
    case MFLOCK_ARGS_OK:
        break;
    case MFLOCK_ARGS_ERR:
        mflock_usage();
        mflock_free(mf);
        return EXIT_FAILURE;
        break;
    case MFLOCK_ARGS_QUIT:
        mflock_free(mf);
        return EXIT_SUCCESS;
        break;
    }

    mflock_init(mf, argc, argv);

    if(mf->verbose > 1)
    {
        mflock_show(mf, stdout);
    }

    mflock_show(mf, mf->logf);

    mflock_run(mf);

    mflock_summary(mf);

    /* Write coordinates to disk */
    mflock_save_coordinates(mf);

    // Write information about which contacts were enabled at the
    // end (when --autopairs is used).
    mflock_write_autopairs(mf);

    /* Write chimera cmm file (.gz) */
    mflock_write_cmm(mf);

    /* Write a few last things and close the log file */
    mflock_close_log(mf);

    /* Free most data */
    mflock_free(mf);

    return EXIT_SUCCESS;
}
