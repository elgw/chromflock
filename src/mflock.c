/**
 * @file mflock.c
 * @author Erik Wernersson
 * @date 2020-2023
 */

#include "mflock.h"

// TODO: Offer alternative for non x86-systems
#include "fast_prng/normal.h"


/* Since mflock is a command line tool it typically fails by just
 * returning with EXIT_FAILURE without freeing all memory.
 *
 * - The error free path should have no memory leaks.
 */

/* Forward declarations */
typedef struct {
    uint32_t * I; // List with pairwise distances
    size_t n_pairs; // Number of pairs in I
    double * beads; /* Bead coordinates */
    size_t n_beads; // number of points
    int diploid; /* Cast the labels to diploid format */
    uint8_t * L; // chr labels per bead
    double * R; // wanted radii together with kRad

    // Geometry
    // Sphere if E isn't set.
    double r0; // bead radius
    double volq;
    elli * E; // Ellipse

    // Initialization
    size_t rseed;

    // For the optimizer
    size_t maxiter;
    size_t maxtime;

    // From the optimizer
    double grad_final;
    size_t iter_final;
    size_t time_final;
    double err_final;

    // Input file names
    char * wfname; // File with contact indications (depreciated)
    char * contact_pairs_file; /* File with contact pairs */
    char * rfname; // Wanted radius (GPSeq)
    char * xfname; // Initial X
    char * lfname; // Chromosome labels

/* one well per bead data */
    char * fname_bead_wells;
    uint32_t n_bead_wells;
    double * bead_wells;

    int newx;

    // Output file names
    char * xoutfname;
    char * ofoldername; // Outfolder
    char * logfname; // log file name
    int cmmz;

    FILE * logf;

    // general
    int verbose;
    double compress; // compress chromosomes by attracting them to their COMs
    int liveView;

    /* Name of lua script to handle the beads dynamics */
    char * luaDynamicsFile;
} mflock_t;




/** @brief create a new default configuration
 * free with mflock_free
*/
static mflock_t * mflock_new(void);

/** @brief free an mflock_t
 *
 * also frees everything that it points to
 */
static void mflock_free(mflock_t * p);

/** @brief Print the settings to FILE
 *
 * f can of course be stdout.
 */
static void mflock_show(mflock_t * p, FILE * f);

/** @brief Report status of mflock to log and screen
 *
 */
static void mflock_summary(mflock_t * p);

/** @brief Read contact pairs from a binary file
    Sets p->I (the contacts) and p->NI (number of contact pairs)
    @return - Nothing, but aborts the program on failure.
*/
static void mflock_read_contact_pairs(mflock_t * p);

/** @brief Load radial constraints
 *
 * only if rfname is set
 * Read GPSeq radius values as binary double.
 * If that does not work, try as text, one value per line
 */
static int mflock_load_radial_constraints(mflock_t * p);

/** @brief Load or set new coordinates
 *
 * Tries to call mflock_init_coordinates
 *
*/

static void mflock_init_coordinates(mflock_t * p);


/** @brief Load bead coordinates from csv file
 *
 */
static int mflock_load_coordinates(mflock_t * p);


/** @brief Write coordinates to disk, also write the column names
 * to the log file  */
static void mflock_save_coordinates(mflock_t * p);

/** @brief Read label matrix pointed to by p->lfname
 */
static int mflock_load_bead_labels(mflock_t * p);

/* For logging */
static void mflock_logwrite(const mflock_t * p, int level, const char * fmt, ...);

/**
 * @breif The beads dynamics main loop
 *
 * @param p the settings
 * @param Fb: the Brownian force
 */
static int mflock_dynamics(mflock_t * restrict p,
                   double Fb);

typedef enum {
    MFLOCK_ARGS_OK,
    MFLOCK_ARGS_ERR,
    MFLOCK_ARGS_QUIT
} mflock_cli_status;

/** @brief parse command line arguments */
static mflock_cli_status mflock_parse_cli(mflock_t * p, int argc, char ** argv);

/** @brief initialization from valid command line arguments
*/
static void mflock_init(mflock_t * p, int argc, char ** argv);


/** @brief set the bead size from the volume quotient */
static void mflock_set_bead_size(mflock_t * p);

/** @brief Per Chr Centre of mass compression
 *
 *Apply a force that attracts each bead to the centre of mass of it's
 * chromosome This slows down the computations quite much because the
 * beads get close to each other so the collision detection gets more to do.
 *
 * TODO: Unnecessary to allocate/free things here and to count the
 * number of beads per chromosome.
 */
static void comforce(mflock_t * restrict p,
                     double * restrict G);

/* END OF FORWARD DECLARATIONS */

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


static double norm3(const double * restrict X)
{
    double n = 0;
    for(size_t kk = 0; kk<3; kk++)
        n+=pow(X[kk], 2);
    return sqrt(n);
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


/* Euclidean distance between two 3D-vectors */
static double eudist3p2(const double * A, const double * B)
{
    return pow(A[0]-B[0], 2) + pow(A[1]-B[1], 2) + pow(A[2]-B[2],2);
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


static int mflock_dynamics(mflock_t * restrict p,
                           double Fb)
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
    if(p->E != NULL)
    {
        /* add ellipse parameters otherwise sphere domain */
        fconf.E = p->E;
        fconf.Es = elli_new(
                            p->E->a - fconf.r0,
                            p->E->b - fconf.r0,
                            p->E->c - fconf.r0);
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
    double Frnd = .10; /* Don't change this, change Fb instead */

    lua_State *L = NULL;


    L = luaL_newstate();
    luaL_openlibs(L); // might not be needed, performance penalty?
    lua_register(L, "usleep", usleep_lua);

    if (luaL_dofile(L, p->luaDynamicsFile))
        luaerror(L, "cannot run config. file: %s", lua_tostring(L, -1));


    int luaquit = 0;
    size_t iter = 0;
    do{
        iter++;
        p->iter_final++;

        /* Get settings from the lua script */

        lua_getglobal(L, "getConfig"); /* function to be called */
        lua_pushnumber(L, iter); /* push arguments */
        lua_pushnumber(L, p->newx);

        /* do the call (2 arguments, 0 result) */
        if (lua_pcall(L, 2, 0, 0) != LUA_OK)
            luaerror(L, "error running function 'getConfig' in %s : %s",
                     p->luaDynamicsFile,
                     lua_tostring(L, -1));

        /* retrieve result */
        fconf.kDom = lua_get_float(L, "kDom");
        fconf.kVol = lua_get_float(L, "kVol");
        fconf.kInt = lua_get_float(L, "kInt");
        fconf.kRad = lua_get_float(L, "kRad");
        fconf.kBeadWell = lua_get_float(L, "kBeadWell");
        fconf.kChrWell = lua_get_float(L, "kChrWell");
        p->compress = lua_get_float(L, "kCom");
        Fb = lua_get_float(L, "fBrown");
        luaquit = lua_get_int(L, "quit");
        fconf.dInteraction = fconf.r0 * lua_get_float(L, "dInteraction");

        if(p->verbose > 10)
        {
            printf("kInteraction: %f\n", fconf.kInt);
            printf("dInteraction: %f\n", fconf.dInteraction);
        }

        /* Start of Molecular Dynamics
         * 2.1 Gradient of functional */
        grad3(X,
              p->n_beads,
              p->R,
              p->I,
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
                double d[] = {0,0,0};
                d[0] = normal(); /* From normal distribution */
                d[1] = normal();
                d[2] = normal();

                for(size_t idx = 0 ; idx<3; idx++)
                {
                    g[3*kk+idx] += Fb*Frnd*d[idx];
                }
            }
        }

        /* Bead wells */
        if( p->n_bead_wells > 0 )
        {
            bead_wells_gradient(&fconf, p->bead_wells, p->n_bead_wells, X, g);
        }

        /* 2.3 Dampening */
        /* Estimate velocities, note: per component */

        for(size_t pp = 0 ; pp < 3*p->n_beads ; pp++)
        {
            v[pp] = (X[pp] - Xm[pp]) / (2.0 * dt);
        }

        for(size_t kk = 0; kk < 3*p->n_beads; kk++)
        {
            g[kk] = g[kk] + damp*v[kk];
        }

        /* 3. Update X */
        for(size_t pp = 0 ; pp < 3*p->n_beads ; pp++)
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
        double pr = norm3(X+3*kk);
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
    mflock_logwrite(p, 2,  "   Max:  % .3f, % .3f, % .3f, % .3f\n", max, may, maz, mar);
    mflock_logwrite(p, 2,  "   Mean: % .3f, % .3f, % .3f, % .3f\n", mex, mey, mez, mer);
    mflock_logwrite(p, 2,  "   Min:  % .3f, % .3f, % .3f, % .3f\n", mix, miy, miz, mir);

    if(run == 0)
    {
        mflock_logwrite(p, 0, "abnormal exit (Ctrl+c was pressed?)\n");
    }

    char * timestr = cf_timestr();
    mflock_logwrite(p, 2, "Finished at: %s\n", timestr);
    free(timestr);
    return;
}


static void mflock_read_contact_pairs(mflock_t * p)
{
    mflock_logwrite(p, 1, "Reading pairwise interactions from %s\n",
             p->contact_pairs_file);
    uint64_t nCP = 0;
    p->I = contact_pairs_read(p->contact_pairs_file, &nCP);
    if(p->I == NULL)
    {
        fprintf(stderr, "%s/%d Failed to read contact pairs from %s\n",
                __FILE__, __LINE__, p->contact_pairs_file);
        exit(EXIT_FAILURE);
    }
    p->n_pairs = nCP;
    mflock_logwrite(p, 1, "Read %lu contacts pairs\n", nCP);
    return;
}

static void mflock_init_coordinates(mflock_t * p)
{
    p->newx = 0;
    if(p->xfname != NULL)
    {
        p->beads = malloc(3*p->n_beads*sizeof(double));
        assert(p->beads != NULL);
        if(mflock_load_coordinates(p) != 0)
        {
            mflock_logwrite(p, 2, "Could not open x-file, starting from random\n");
            free(p->beads);
            p->beads = NULL;
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
                if(elli_getScale(E, p->beads+3*kk) <= 1)
                { accepted = 1;}
            }
        }

        if(p->E == NULL)
        {
            free(E);
        }
        mflock_logwrite(p, 1, "X[0] = %f\n", p->beads[0]);
    }
    return;
}

static void mflock_validate_labels(const mflock_t * p)
{
    int ok = 1;
    for(size_t kk = 0; kk< p->n_beads; kk++)
    {
        if(p->L[kk] > 31 || p->L[kk] == 0)
        {
            ok = 0;
        }
    }
    if( ! ok )
    {
        mflock_logwrite(p, 0,
                 "The labels are not as expected. Any label l "
                 "should satisfy 0 < L < 32\n");
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


static void mflock_save_coordinates(mflock_t * p)
{
    const double * restrict X = p->beads;
    size_t N = p->n_beads;

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
        mflock_logwrite(p, 1, "Columns: x, y, z, r");

        mflock_logwrite(p, 1, "\n");
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

        fprintf(f, "\n");
    }
    fclose(f);
    return;
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
    fprintf(f, "    - Model parameters\n");
    fprintf(f, "    Dynamics program: %s\n", p->luaDynamicsFile);
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

    fprintf(f, "    diploid=%d\n", p->diploid);

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

static void mflock_usage()
{
    printf("mflock %s usage:\n"
           "\n", cf_version);
    printf("Required arguments:\n");
    printf(" --contact-pairs <file>\n\t"
           "File with contact pairs (uint32).\n");
    printf(" --labels <file>, --L <file>\n\t"
           "bead label array (uint8_t). Determines the number of beads.\n\t"
           "(2X if --diploid is set)\n");
    printf(" --config mflock.lua\n\t"
           "lua script that configures the dynamics dynamics\n");
    printf("\n"
           "Optional arguments:\n");
    printf(" --coordinates <file>, -x <file>\n\t"
           "file with initial coordinates\n");
    printf(" --diploid\n\t"
           "Interpret/cast the labels for a diploid structure by duplication\n");
    printf(" --seed N, -s N\n\t"
           "seed for random number generator (defaults to random)\n");
    printf(" --radii <file>, -r <file>\n\t"
           "file with wanted radii \n"
           "\tall files should be encoded as 64-bit floats\n");
    printf(" --maxiter N, -n N\n\t"
           "maximum number of iterations\n");
    printf(" --maxtime N, -t N\n\t"
           "time budget in seconds\n");
    printf(" --bead-wells <file>\n\t"
           "specify a file containing preferred bead locations (bead wells)\n");
    printf(" --config-show mflock.lua\n\t"
           "Perform a dry run and print the dynamics settings as a table.\n\t"
           "No additional arguments are required or used\n");
    printf("\n"
           "Geometry\n");
    printf(" --radius r0, -R r0\n\t"
           "bead radius (sphere has radius 1)\n");
    printf(" --volq vq, -Q vq\n\t"
           "volume quotient (bead volume/domain volume).\n\t"
           "For ellipsoidal geometry, set 1=ea>=eb>=ec>0\n");
    printf(" --ea ea\n\t"
           "Major axis of ellipsoid\n");
    printf(" --eb eb\n\t"
           "Second axis of ellipsoid\n");
    printf(" --ec ec\n\t"
           "Third axis of ellipsoid\n");
    printf("\n"
           "General settings\n");
    printf(" --verbose level, -v level\n\t"
           "verbosity (lowest=0, default=1, ... )\n");
    printf(" --outFolder name, -o name\n\t"
           "where to store results\n");
    printf(" --cmmz, -z\n\t"
           "write compressed cmm files (gzip)\n");
    printf(" --live\n\t"
           "enable live monitoring (only when compiled with SDL)\n");
    printf(" --help, -h\n\t"
           "show this help message. For more info see 'man mflock'\n");
    printf(" --defaults, -d\n\t"
           "show default settings for the parameters\n");
    printf("\n");
    printf("For additional help see the man page or visit\n"
        "https://www.github.com/elgw/chromflock/\n");
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
        /* Data */
        { "wFile",         required_argument, NULL,   'w' },
        { "contact-pairs", required_argument, NULL,   'p' },
        { "xFile",         required_argument, NULL,   'x' },
        { "coordinates",   required_argument, NULL,   'x' },
        { "rFile",         required_argument, NULL,   'r' },
        { "radii",         required_argument, NULL,   'r' },
        { "labels",        required_argument, NULL,   'L' },
        { "lFile",         required_argument, NULL,   'L' },
        { "outFolder",     required_argument, NULL,   'o' },
        { "bead-wells",    required_argument, NULL,   'W' },
        /* Settings */
        { "diploid",       no_argument,       NULL,   'D' },
        { "maxiter",       required_argument, NULL,   'n' },
        { "maxtime",       required_argument, NULL,   't' },
        { "seed",          required_argument, NULL,   's' },
        { "verbose",       required_argument, NULL,   'v' },
        { "live",          no_argument,       NULL,   'a' },
        { "cmmz",          no_argument,       NULL,   'z' },
        { "defaults",      no_argument,       NULL,   'd' },
        /* Geometry */
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
        { NULL,            0,                 NULL,    0  }
    };

    int ch;
    while((ch = getopt_long(argc, argv,
                            "A:B:C:Dw:x:r:n:t:R:v:o:p:hMs:L:zcadQ:l:W:",
                            longopts, NULL)) != -1)
    {
        switch(ch) {
        case 'i':
            printf("mflock (chromflock version %s)\n", cf_version);
            printf("Build date: %s, %s\n", __DATE__, __TIME__);
            printf("GIT HASH: %s\n", GIT_VERSION);
            printf("Compiler: %s\n", CC_VERSION);
            return MFLOCK_ARGS_QUIT;
        case 'd':
            printf("Defaults:\n");
            mflock_show(p, stdout);
            return MFLOCK_ARGS_QUIT;
        case 'D':
            p->diploid = 1;
            break;
        case 'w':
            free(p->wfname);
            p->wfname = strdup(optarg);
            assert(p->wfname != NULL);
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
            assert(p->contact_pairs_file != NULL);
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
        case 'z':
            p->cmmz = 1;
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
            return MFLOCK_ARGS_ERR;
        }
    }

    if(p->luaDynamicsFile == NULL)
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
        if(p->ofoldername[strlen(p->ofoldername) - 1] != '/' )
        {
            p->ofoldername = realloc(p->ofoldername, strlen(p->ofoldername)+1);
            p->ofoldername[strlen(p->ofoldername)+1] = '\0';
            p->ofoldername[strlen(p->ofoldername)] = '/';
        }
    }

    if(p->contact_pairs_file == NULL)
    {
        fprintf(stderr, "WARNING: Contacts file not set\n");
        return MFLOCK_ARGS_ERR;
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

static mflock_t *  mflock_new(void)
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
        size_t bead1 = mf->bead_wells[4*kk]+1;
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

static void mflock_init(mflock_t * mf, int argc, char ** argv)
{
    mflock_set_and_create_output_folder(mf);

    /* Set names of output files */
    mf->xoutfname = malloc(1024*sizeof(char));
    assert(mf->xoutfname != NULL);
    sprintf(mf->xoutfname, "%s%s", mf->ofoldername, "coords.csv");

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

    /* Once we know how many beads we can set their radius based on the
     * volume quotient (or do nothing if it was given at the command line) */
    mflock_set_bead_size(mf);

    /* Read coordinates from previous iteration, or generate
       random coordinates if first iteration.
    */
    mflock_init_coordinates(mf);

    mflock_read_contact_pairs(mf);

    mflock_load_radial_constraints(mf);

    mflock_load_bead_wells(mf);

    return;
}


void mflock_free(mflock_t * p)
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
    free(p->beads);
    free(p->fname_bead_wells);
    free(p->bead_wells);
    free(p);
    return;
}


static void mflock_set_bead_size(mflock_t * p)
{
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
        p->r0 = cbrt( 3.0*vq*Vd / (4.0*p->n_beads*M_PI) );

        if(p->verbose > 2)
        {
            double Vbeads = p->n_beads*pow(p->r0,3)*M_PI*4.0/3.0;
            printf("Vd: %f, Vbeads: %f, Vbeads/Vd = %f\n", Vd, Vbeads, Vbeads/Vd);
            assert(fabs(Vbeads/Vd - vq)<1e-5);
        }
    }
}

/** @brief write to both stdout and to the log file if p->verbose >= level */
static void mflock_logwrite(const mflock_t * p, int level, const char *fmt, ...)
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


static int mflock_load_coordinates(mflock_t * p)
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
    for(size_t ll = 0; ll<p->n_beads; ll++)
    {
        int read = getline(&line, &len, f);
        if(read == -1)
        {
            printf("Failed to read line %zu\n", ll+1);
            return -1;
        }
        char *ptr = strtok(line, delim);
        p->beads[3*ll] = atof(ptr);
        ptr = strtok(NULL, delim);
        p->beads[3*ll+1] = atof(ptr);
        ptr = strtok(NULL, delim);
        p->beads[3*ll+2] = atof(ptr);
    }

    free(line);
    return 0;
}

static void * solve_t(void * args)
{
    mflock_dynamics((mflock_t *) args, 0);
    return NULL;
}

/* Write chimera file (for simple visualization) */
static void mflock_write_cmm(const mflock_t * p)
{
    const double * restrict X = p->beads;
    if(p->cmmz == 1)
    {
        char * cmmfile = malloc(1024*sizeof(char));
        assert(cmmfile != NULL);
        sprintf(cmmfile, "%s/cmmdump.cmm.gz", p->ofoldername);

        cmmwritez(cmmfile, X, p->n_beads, p->r0, p->I, p->n_pairs, p->L);
        free(cmmfile);
    } else {
        char * cmmfile = malloc(1024*sizeof(char));
        assert(cmmfile != NULL);
        sprintf(cmmfile, "%s/cmmdump.cmm", p->ofoldername);
        cmmwrite(cmmfile, X, p->n_beads, p->r0, p->I, p->n_pairs, p->L);
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
        int quit = 0;

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

        liveview(p->beads, p->L, p->n_beads, &quit, p->r0, E);
        pthread_join(th, NULL);
    } else {
        mflock_dynamics(p, 0);
    }
#else
    mflock_dynamics(p, 0);
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

    if(mf->verbose>0)
    {
        mflock_show(mf, stdout);
    }

    mflock_show(mf, mf->logf);

    mflock_run(mf);

    mflock_summary(mf);

    /* Write coordinates to disk */
    mflock_save_coordinates(mf);

    /* Write chimera cmm file (.gz) */
    mflock_write_cmm(mf);

    /* Write a few last things and close the log file */
    mflock_close_log(mf);

    /* Free most data */
    mflock_free(mf);

    return EXIT_SUCCESS;
}
