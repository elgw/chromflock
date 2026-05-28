#include "ellipsoid.h"

typedef double (*elli_distance) (const elli * E,
                                 const double * restrict p,
                                 double * restrict X);

typedef struct {
    const char * name;
    elli_distance fun;
    size_t maxiter;
    double eps;
    double time_s;
} elli_d_method;

static double
eudist3(const double * a, const double * b)
{
    return sqrt(
        pow(a[0]-b[0], 2) +
        pow(a[1]-b[1], 2) +
        pow(a[2]-b[2], 2));
}


static double
norm3(const double * a)
{
    return sqrt(
        pow(a[0], 2) +
        pow(a[1], 2) +
        pow(a[2], 2));
}

static void v3_print(const double * x)
{
    printf("(%f, %f, %f)\n", x[0], x[1], x[2]);
    return;
}

static double sprod(double * a, double * b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}


static double compare_sphere(double r)
{
    double p[] = {.99, .99, .99};
    // Compare with sphere
    elli * E = elli_new(r, r, r);
    assert(fabs(elli_vol(E) - pow(r,3)*4.0/3.0*M_PI) < 1e-9);
    for(size_t kk = 0; kk<10000; kk++)
    {
        p[0] = r*(2.0*(rand()/ (double) RAND_MAX) - 1.0);
        p[1] = r*(2.0*(rand()/ (double) RAND_MAX) - 1.0);
        p[2] = r*(2.0*(rand()/ (double) RAND_MAX) - 1.0);

        // Classifies inside/outside correctly
        double d = norm3(p);
        if(d<r)
            assert(elli_isInside(E, p) == 1);
        if(d>r)
            assert(elli_isInside(E, p) == 0);

        // Distances to surface
        double de = elli_dist_p(E, p, NULL);
        double dg = elli_gdist(E, p, NULL);

        if( fabs(fabs(de) - fabs(r-norm3(p))) > 1e-6)
        {
            printf("elli_dist_p = %f, |1-||p||| = %f\n", de, fabs(r-norm3(p)));
            assert(0);
        }

        if( fabs(fabs(dg) - fabs(r-norm3(p))) > 1e-5)
        {
            printf("elli_gdist = %f, |1-||p||| = %f\n", dg, fabs(r-norm3(p)));
            assert(0);
        }

        // Radius
        double re = elli_radius(E, p);
        if( fabs(re - norm3(p)/r) > 1e-6)
        {
            printf("Radius mismatch:\n");
            elli_print(stdout, E);
            printf("p = %f %f %f\n", p[0], p[1], p[2]);
            printf("re=%f, norm(p)/r=%f\n", re, norm3(p)/r);
            assert(0);
        }

    }
    free(E);
    return 0;
}

void geodesic_random()
{
    elli * E = elli_new(2, 1.1, 1);
    double p[3];
    double r = 5; // max radius
    printf("\n");
    for(size_t kk = 0; kk<10000; kk++)
    {
        //    printf("\r %zu", kk); fflush(stdout);
        p[0] = r*(2.0*(rand()/ (double) RAND_MAX) - 1.0);
        p[1] = r*(2.0*(rand()/ (double) RAND_MAX) - 1.0);
        p[2] = r*(2.0*(rand()/ (double) RAND_MAX) - 1.0);
        elli_gdist(E, p, NULL);
    }
    free(E);
    return;
}

double * getPointsOnSphere(size_t N)
{
    double * P = malloc(3*N*sizeof(double));
    assert(P != NULL);
    for(size_t kk = 0; kk<N; kk++)
    {
        int gotPoint = 0;
        while(gotPoint == 0)
        {
            double * X = P+3*kk;
            for(int ll = 0; ll<3; ll++)
                X[ll] = (2.0*(double) rand()/ (double) RAND_MAX) - 1.0;

            double n = norm3(X);
            if(n<=1 && n>0)
            {
                for(int ll = 0; ll<3 ; ll++)
                    X[ll] /= n;
                gotPoint = 1;
            }
        }
    }
    return P;
}

static double clockdiff(struct timespec* start, struct timespec * finish)
{

    double elapsed = (finish->tv_sec - start->tv_sec);
    elapsed += (finish->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}

void smalltest()
{
    /* Test smaller ellipse as primary distance test
     */

    printf("--> Smaller ellipsoid as lower primary distance check\n");

    srand(time(NULL));
    double delta = 0.02;
    double a = 1.7; double b = 1; double c = 1;
    elli * E = elli_new(a, b, c);
    elli * Ei = elli_new(a - delta, b - delta, c - delta);
    printf("    E = %f %f %f\n", E->a, E->b, E->c);
    printf("    Ed = %f %f %f  (E - %f)\n", Ei->a, Ei->b, Ei->c, delta);
    printf("    Testing that f(p) < d for p in Ed\n");

    E->maxiter = 1000;
    for(size_t kk = 0; kk<20000; kk++)
    {
        double X[3];
        double Y[3];
        X[0] = (double) rand()/ (double) RAND_MAX;
        X[1] = (double) rand()/ (double) RAND_MAX;
        X[2] = (double) rand()/ (double) RAND_MAX;
        elli_project(Ei, X); // scale X to be on Ei
        double d = elli_gdistR(E, X, Y);
        if(d>=delta)
        {
            printf("No, for point %zu:", kk);
            printf("delta: %f d(xi, E) = %f\n", delta, eudist3(X,Y));
            assert(0);
        }
    }
    free(E);
    free(Ei);
    return;
}

static void
test_distance_to_ellipsoid(int argc, char ** argv)
{
    assert(argc == 7);

    double X[6];
    for(int kk = 1; kk<7; kk++)
    {
        X[kk-1] = atof(argv[kk]);
    }
    elli * E = elli_new(X[0], X[1], X[2]);
    double Y[3];
    double d = elli_gdist(E, X+3, Y);
    printf("d = %.10f\n", d);
    printf("Y = [%.10f %.10f %.10f]\n", Y[0], Y[1], Y[2]);
}


// Test speed vs precision for a specific ellipse
// N points on the ellipsoid are generated. Then they are slightly
// displaced along the normal and it is checked how well the
// methods find back the original point.
static void
test_speed_vs_precision(int argc, char ** argv)
{
    assert(argc == 4);
    size_t N = 1e6;
    elli * E = elli_new(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));

    double dist = 0.02*E->c; // distance from E

    E->eps = 1e-5;
    printf("eps: %e\n", E->eps);

    printf("Generating %zu random points on a sphere\n", N);
    double * P = getPointsOnSphere(N);
    assert(P != NULL);

    printf("Moving points from sphere to ellipse\n");
    for(size_t kk = 0; kk<N; kk++){
        elli_project(E, P + 3*kk);}

    printf("Generating (exterior) reference points\n");
    double * Q = malloc(3*N*sizeof(double));
    assert(Q != NULL);
    for(size_t kk = 0; kk<N; kk++) {
        double normal[3];
        elli_nnormal(E, P+3*kk, normal);
        for(size_t ll = 0; ll<3; ll++) {
            Q[3*kk+ll] = P[3*kk+ll] + dist*normal[ll];
        }
        assert(elli_isInside(E, Q+3*kk) == 0);
    }

    // Output points
    double * Y = malloc(3*N*sizeof(double));
    assert(Y != NULL);

    elli_d_method methods[] = {
        {"projection", elli_dist_p, 0, 0, 0},
        {"    Bektas", elli_gdist,  2, 1e-9, 0},
        {"      Hart", elli_gdistR, 2, 1e-9, 0},
        {"  Lagrange", elli_gdistL, 2, 1e-9, 0},
        {"    Bektas", elli_gdist,  3, 1e-9, 0},
        {"      Hart", elli_gdistR, 3, 1e-9, 0},
        {"  Lagrange", elli_gdistL, 3, 1e-9, 0},
        {"    Bektas", elli_gdist,  9, 1e-9, 0},
        {"      Hart", elli_gdistR, 9, 1e-9, 0},
        {"  Lagrange", elli_gdistL, 9, 1e-9, 0},
        {NULL, NULL, 0, 0, 0}
    };

    printf("%12s, %10s, %10s, %10s, %10s,  %10s\n",
           "Method",
           "max(abs(err))",
           "mean(abs(err))",
           "tangerr",
           "time [s]",
           "points/s");

    int method_id = 0;
    while(methods[method_id].fun != NULL)
    {
        elli_d_method method = methods[method_id];
        elli_distance dfun = method.fun;
        method_id++;

        E->maxiter = method.maxiter;
        E->eps = method.eps;

        struct timespec tstart;
        struct timespec tend;
        clock_gettime(CLOCK_REALTIME, &tstart);
        for(size_t kk = 0; kk<N; kk++)
        {
            double d = dfun(E, Q+3*kk, Y+3*kk);
            //assert(d > 0);
        }
        clock_gettime(CLOCK_REALTIME, &tend);

        double maxerror = 0;
        double maxnerror = 0;
        double mean_error = 0;
        for(size_t kk = 0; kk<N; kk++)
        {
            double error = eudist3(P+3*kk, Y+3*kk);
            mean_error += fabs(error);
            double nP[3];
            double nPE[3];
            elli_nnormal(E, P+3*kk, nP);
            elli_nnormal(E, Y+3*kk, nPE);
            double nerror = 360.0/M_PI*acos(sprod(nP, nPE));
            if(nerror > maxnerror)
                maxnerror = nerror;

            if(error > maxerror)
                maxerror = error;
        }
        mean_error /= (double) N;
        double ttotal = clockdiff(&tstart, &tend);

        printf("%s-%d,      %.2e,       %.2e,   %.2e,   %.2e,    %.2e\n",
               method.name, (int) method.maxiter,
               maxerror, mean_error,
               maxnerror,
               ttotal, (double) N / ttotal);
    }

    free(Q);
    free(E);
    free(Y);
    return;
}

void v3_random(double X[3], double lim)
{
    for(int kk = 0; kk < 3; kk++)
    {
        X[kk] = lim* (2.0 * (((double) rand() / (double) RAND_MAX)) - 0.5);
    }
}

void
test_gdistS()
{
    elli * E = elli_new(1, 1, 1);
    elli_print(stdout, E);
    double P[3];
    double Y[3] = {0, 0, 0};
    v3_random(P, 1);
    double d = elli_gdistS(E, P);
    double dref = elli_gdist(E, P, Y);
    printf("P = (%f, %f, %f)\n", P[0], P[1], P[2]);
    printf("inside: %d\n", elli_isInside(E, P));
    printf("||P|| = %f\n", norm3(P));
    printf("gdistS: %f\n", d);
    printf("gdist: %f\n", dref);



    free(E);
}

void
test_elli_normal(void) {
    printf("-> elli_normal\n");
    elli * E = elli_new(1, 1, 1);
    double n[3];
    double p[3] = {1, 0, 0};
    elli_normal(E, p, n);
    assert(n[0] != 0.0);
    assert(n[1] == 0.0);
    assert(n[2] == 0.0);
    free(E);
}

void
test_elli_nnormal(void) {
    printf("-> elli_nnormal\n");
    elli * E = elli_new(2, 2, 2);
    double n[3];
    double p[3] = {2, 0, 0};
    double r[3] = {1, 0, 0};
    elli_normal(E, p, n);
    if(eudist3(r, n) > 1e-6 )
    {
        printf("  Normal="); v3_print(n);
        printf("Expected="); v3_print(r);
        assert("Wrong normal calculated");
    }
    free(E);
}

int main(int argc, char ** argv)
{
#ifdef NDEBUG
    fprintf(stderr,
            "Please turn on debugging or some checks will not be performed!\n"
        );
#endif

    if(argc > 1) {
        if(strcmp(argv[1], "--help") == 0) {
            printf("Usage:\n");
            printf("%s"
                   "\n\t run some basic tests and quit"
                   "\n",
                   argv[0]);
            printf("%s A B C"
                   "\n\tTest with an ellipsoid with major axes A, B, C"
                   "\n\tand compare the geodesic distance implementations"
                   "\n",
                   argv[0]);
            printf("%s A B C x y z"
                   "\n\tCompute the distance from (x,y,z) to the ellipoid defined by (A, B, C)"
                   "\n",
                   argv[0]);
            return EXIT_SUCCESS;
        }
    }


    // test_gdistS(); -- test this with the other ones

    test_elli_normal();
    test_elli_nnormal();

    if(argc == 4)
    {
        test_speed_vs_precision(argc, argv);
        exit(EXIT_SUCCESS);
    }

    if(argc == 7)
    {
        test_distance_to_ellipsoid(argc, argv);
        exit(EXIT_SUCCESS);
    }

    smalltest();

    for(double r = .8; r<1.2; r=r+0.01)
    {
        compare_sphere(r);
    }

    geodesic_random();

    printf("All tests passed\n");
    return 0;
}
