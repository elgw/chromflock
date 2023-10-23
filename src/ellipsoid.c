#include "ellipsoid.h"

static double eudist3(const double * a, const double * b)
{
    return sqrt(
        pow(a[0]-b[0], 2) +
        pow(a[1]-b[1], 2) +
        pow(a[2]-b[2], 2));
}

/*
  static void  cross3(double * a, double * b, double * C)
  {
  C[0] = (a[1]*b[2]-a[2]*b[1]);
  C[1] = - (a[0]*b[2]-a[0]*b[2]);
  C[2] = (a[0]*b[1]-a[1]*b[0]);
  return;
  }
*/

static double norm3(const double * a)
{
    return sqrt(
        pow(a[0],2) +
        pow(a[1],2) +
        pow(a[2],2));
}

elli * elli_new(double a, double b, double c)
{
    elli * E = malloc(sizeof(elli));
    assert(E != NULL);
    E->a=a;
    E->a2=pow(a,2);
    E->b=b;
    E->b2=pow(b,2);
    E->c=c;
    E->c2=pow(c,2);
    E->maxiter = 40;
    E->eps = 1e-9;
    return E;
}

void elli_show(elli * E)
{
    printf("Ellipse: x^2/a^2 + y^2/b^2 + z^2/c^2 = 1, maxiter: %zu, eps: %e\n", E->maxiter, E->eps);
    printf("         a=%f, b=%f, c=%f, volume=%f\n",
           E->a, E->b, E->c, elli_vol(E));
    return;
}

int elli_isInside(elli * E, double * p)
{
    double t =
        pow(p[0],2)/E->a2 +
        pow(p[1],2)/E->b2 +
        pow(p[2],2)/E->c2;

    if(t<=1)
        return 1;

    return 0;

}

double elli_vol(elli * E)
{
    return 4.0/3.0*M_PI*E->a*E->b*E->c;
}

double elli_radius(elli * E, double * p)
{

    double t = sqrt( 1.0/(  pow(p[0]/E->a, 2)
                            + pow(p[1]/E->b, 2)
                            + pow(p[2]/E->c, 2)
                         ));

    return 1.0/t;
}

double elli_pdist(elli * E, double * p)
{

    // Distance from point to surface
    // by projection through origo
    double t = sqrt( 1.0/(  pow(p[0]/E->a, 2)
                            + pow(p[1]/E->b, 2)
                            + pow(p[2]/E->c, 2)
                         ));

    double ip[] = {0,0,0};
    for(int kk = 0; kk<3; kk++)
        ip[kk] = t*p[kk];

    return eudist3(p, ip);
}

double elli_getScale(const elli * restrict E, const double * restrict X)
{
    return sqrt(elli_getScale2(E, X));
}

double elli_getScale2(const elli * restrict E, const double * restrict X)
{
    return pow(X[0]/E->a,2) + pow(X[1]/E->b,2) + pow(X[2]/E->c,2);
}

void elli_scale(elli * E, double * X)
{
    double z = sqrt(pow(X[0]/E->a,2) + pow(X[1]/E->b,2) + pow(X[2]/E->c,2));
    if(z == 0)
    {
        X[0] = 0; X[1] = 0; X[2] = 0;
    } else {
        X[0]*=1/z; X[1]*=1/z; X[2]*=1/z;
    }
}

void elli_scale2(elli * E, double * X)
{
    // Two-step scaling
    double X0[3];
    memcpy(X0, X, 3*sizeof(double));
    elli_scale(E, X);
    double n[3];
    elli_nnormal(E, X, n);
    double d = eudist3(X0, X);
    double X1[3];
    for(int ll = 0; ll<3; ll++)
        X1[ll] = X[ll] + d*n[ll];
    for(int ll = 0; ll<3; ll++)
        X[ll] += (X0[ll]-X1[ll]);

}

void elli_nnormal(elli * E, double * X, double * N)
{
    N[0] = 2.0*X[0]/pow(E->a,2);
    N[1] = 2.0*X[1]/pow(E->b,2);
    N[2] = 2.0*X[2]/pow(E->c,2);
    double n = norm3(N);
    for(int kk = 0; kk<3; kk++)
    {
        N[kk]/=n;
    }
}

void elli_normal(elli * E, double * X, double * N)
{
    N[0] = 2.0*X[0]/pow(E->a,2);
    N[1] = 2.0*X[1]/pow(E->b,2);
    N[2] = 2.0*X[2]/pow(E->c,2);
}

static double f_div_df(const double l,
                       const double * restrict E,
                       const double * restrict P)
{
    double f = 0;
    double df = 0;

    for(int kk = 0; kk<3; kk++)
    {
        double pkk2 = pow(P[kk],2);
        double lekk12 = pow(l*E[kk] + 1, 2);
        f+=E[kk]*pkk2/lekk12;
        df+= pow(E[kk],2)*pkk2/pow(l*E[kk] + 1.0, 1)/lekk12;
    }

    df= - 2.0*df;
    f = f - 1.0;

    return f/df;
}

static double elli_dgistR_fun(double * E, double * P, double alpha)
{
    double a = E[0];
    double b = E[1];
    double c = E[2];
    double x2 = pow(P[0], 2);
    double y2 = pow(P[1], 2);
    double z2 = pow(P[2], 2);
    double pa = pow(a + alpha/a, 2);
    double dpa  = 2.0*alpha/pow(E[0],2)+2;
    double pb = pow(b + alpha/b, 2);
    double dpb  = 2.0*alpha/pow(E[1],2)+2;
    double pc = pow(c + alpha/c, 2);
    double dpc  = 2.0*alpha/pow(E[2],2)+2;
    double f = pa*pb*pc - (pb*pc*x2 + pa*pc*y2 + pa*pb*z2);
    double df = (dpa*pb*pc + pa*dpb*pc + pa*pb*dpc) -
        ((dpb*pc + pb*dpc)*x2 +
         (dpa*pc + pa*dpc)*y2 +
         (dpa*pb + pa*dpb)*z2);
    return f/df;
}

double elli_gdistR(elli * restrict ellipse, double * restrict P, double * restrict X)
{
    /* Root finding method, see John C. Hart, 1994 */
    double E[3];
    E[0] = ellipse->a;
    E[1] = ellipse->b;
    E[2] = ellipse->c;

    size_t maxiter = ellipse->maxiter;
    double eps = ellipse->eps;
    size_t iter = 0;
    // Graphics Gems 4
    double alpha = 0;
    double alpha0 = norm3(P)*E[0]; // a>=b>=c
    while(fabs(alpha-alpha0)>eps && iter<maxiter)
    {
        alpha0 = alpha;
        alpha = alpha0 - elli_dgistR_fun(E, P, alpha0);
        iter++;
    }
//printf("%f\n", fabs(alpha-alpha0));
    double a2 = ellipse->a2;
    double b2 = ellipse->b2;
    double c2 = ellipse->c2;

    X[0] = a2*P[0]/(alpha + a2);
    X[1] = b2*P[1]/(alpha + b2);
    X[2] = c2*P[2]/(alpha + c2);

    return 0;
}

double elli_gdistL(elli * restrict ellipse,
                   const double * restrict P,
                   double * restrict Y)
{

    const double eps = ellipse->eps;
    double l0 = 0;
    double l1 = 5e-2; // check typical lambda if solving similar problems over and over
    // Possibly re-use last lambda.
    size_t iter = 0;
    const size_t maxiter = ellipse->maxiter;
    double E[3];
    E[0] = 1.0/ellipse->a2;
    E[1] = 1.0/ellipse->b2;
    E[2] = 1.0/ellipse->c2;
    assert(fabs(ellipse->c2 - pow(ellipse->c,2))<1e-9);

    /* Improved starting condition */
//  l1 = 1.0;

    while( (fabs(l1-l0) > eps) && (iter < maxiter) )
    {
        iter++;
        l0 = l1;
        l1 = l0 - f_div_df(l0, E, P);
    }

    for(size_t kk = 0 ; kk<3; kk++)
    {
        Y[kk] = P[kk]/(1.0 + l1*E[kk]);
    }

//  printf("%e\n", l1);
    if(norm3(P)>norm3(Y))
    {
        return eudist3(P, Y);
    } else {
        return -eudist3(P,Y);
    }
}

double elli_gdist(elli * ellipse, double * P, double * Y)
{
    /* Following Bektas 2015
     * Using The Newton-Raphson method.
     *
     * Note that for interior points, there are potentially infinite number of solutions
     * due to the way that the problem is posed. I.e., only looking at the direction of the normals.
     * The solution is only guaranteed to be a shortest distance if P is outside of the ellipse
     *
     * See also this library:
     * https://tcg.mae.cornell.edu/pubs/Pope_FDA_08.pdf
     */

    // Most iterations to be used
    const size_t maxiter = ellipse->maxiter; // 40
    // Abort when max(fabs(deltaX)) is less than this value
    const double tol = ellipse->eps;

    double a = ellipse->a;
    double b = ellipse->b;
    double c = ellipse->c;

    assert(a>=b); assert(b>=c);

    double E = copysign(1.0, ellipse->a)/ellipse->a2;
    double F = copysign(1.0, ellipse->b)/ellipse->b2;
    double G = copysign(1.0, ellipse->c)/ellipse->c2;

    double X[3]; // Solution at current iteration
    double D[3];

    /*
     *        Initialization
     */

    const char initialization = 'C';
    if(initialization == 'A') {
        /* As in the paper -- fastest, but would fail to converge more often than the other alternatives
         * without the projection down to the ellipsoid */
        double normP = norm3(P);
        X[0] = P[0]*a/normP;
        X[1] = P[1]*b/normP;
        X[2] = P[2]*c/normP;
    }
    if(initialization == 'B') {
        /* This is what was used in the MATLAB code. Possibly because the 'A' was
         * more prone to cause the convergence to fail */
        X[0] = P[0];
        X[1] = P[1];
        X[2] = P[2];
    }
    if(initialization == 'C') {
        /* Project directly on the ellipsoid -- requires more iterations than B */
        double zP = sqrt(pow(P[0]/a,2) + pow(P[1]/b,2) + pow(P[2]/c,2));
        X[0]=P[0]/zP; X[1]=P[1]/zP; X[2]=P[2]/zP;
    }


    double maxdelta = 0; // book keeping
    size_t iter=0;

    for(iter = 0 ; iter<maxiter; iter++)
    {

        /* 7 of the 9 elements of A are non-zeros */
        double Aa = F*X[1] - (X[1]-P[1])*E;
        double Ab = (X[0]-P[0])*F - E*X[0];
        double Ac = G*X[2] - (X[2]-P[2])*E;
        double Ad = (X[0]-P[0])*G - E*X[0];
        double Ae = 2*E*X[0];
        double Af = 2*F*X[1];
        double Ag = 2*G*X[2];

        if(0){
            printf("A = [[%f %f %f];[%f %f %f];[%f %f %f]]\n", Aa, Ab, 0.0, Ac, 0.0, Ad, Ae, Af, Ag);
        }

        // Ab in the .m
        D[0] = (X[0]-P[0])*F*X[1] - (X[1]-P[1])*E*X[0];
        D[1] = (X[0]-P[0])*G*X[2] - (X[2]-P[2])*E*X[0];
        D[2] = E*pow(X[0], 2) + F*pow(X[1], 2) + G*pow(X[2], 2) - 1;

        //   printf("D = [%f %f %f]\n", D[0], D[1], D[2]);

        /*  Wolfram Alpha is my friend
         * inv [[a,b,0],[c,0,d],[e,f,g]] =
         * 1/(a d f + b c g - b d e)[[d f, b g, -b d],
         * [c g - d e,  -a g, a d],
         * [-c f, a f - b e, b c)]]
         */
        double dx =           Ad*Af*D[0]           + Ab*Ag*D[1] - Ab*Ad*D[2];
        double dy = (Ac*Ag - Ad*Ae)*D[0]          -(Aa*Ag)*D[1] + Aa*Ad*D[2];
        double dz =        -(Ac*Af)*D[0] + (Aa*Af - Ab*Ae)*D[1] + Ab*Ac*D[2];
        double f = 1.0/(Aa*Ad*Af + Ab*Ac*Ag - Ab*Ad*Ae);

        dx*=f;    dy*=f;    dz*=f;

        if(0){
            printf("A^-1\n");
            printf("%f , %f, %f\n",           f*Ad*Af,           + f*Ab*Ag, - f*Ab*Ad);
            printf("%f , %f, %f\n", f*(Ac*Ag - Ad*Ae),          -f*(Aa*Ag), + f*Aa*Ad);
            printf("%f , %f, %f\n",        -f*(Ac*Af), + f*(Aa*Af - Ab*Ae), + f*Ab*Ac);
            printf("delta = [%f %f %f]\n", dx, dy, dz);
        }

        // Update X
        X[0] -= dx; X[1] -= dy; X[2] -= dz;

        /* Project onto the ellipsoid */
        // Fixes stability issues for some starting points
//    elli_scale(ellipse, X);


        if(0){
            printf("X = [%f, %f, %f]\n", X[0], X[1], X[2]);
        }

        maxdelta = fabs(dx);
        if(fabs(dy)>maxdelta)
            maxdelta = fabs(dy);
        if(fabs(dz)>maxdelta)
            maxdelta = fabs(dz);
        if(maxdelta < tol)
            iter = maxiter;

        //    printf("Iter: %zu Delta: %f (%f, %f, %f)\n", iter, maxdelta, X[0], X[1], X[2]);
    }

    double dist = copysign(1.0, X[2]-P[2])*copysign(1.0, X[2])*eudist3(P, X);

#ifndef NDEBUG
    /* Verify that the result is correct by predicting
     * the position of the input point based on X and dist:
     *
     * Q: = X + dist*N == P ?
     *
     * N is the surface normal at X.
     * NOTE: ||Q-P|| could be used as a convergence criterion.
     */

    double N[3];
    elli_nnormal(ellipse, X, N);


    double Q[3];
    double normN = norm3(N);
    for(int kk = 0; kk<3; kk++) {
        Q[kk] = X[kk] - dist*N[kk]/normN;
    }

    double error = eudist3(P,Q);

    if(error > 1e-2 && 0)
    {
        printf("\n");
        printf("iter: %zu\n", iter);
        printf("maxdelta: %f\n", maxdelta);
        printf("x^2/a^2 + ... = %f \n", pow(X[0]/a,2) + pow(X[1]/b,2) + pow(X[2]/c,2));
        printf("Error: %e \n", error);
        printf("Ellipse: %f %f %f\n", a, b, c);
        printf("P = [%.10f %.10f %.10f]\n norm(P) = %f\n", P[0], P[1], P[2], norm3(P));
        printf("Q = [%f %f %f]\n norm(Q) = %f\n", Q[0], Q[1], Q[2], norm3(Q));
        printf("%f %f %f %f %f %f\n", a, b, c, E, F, G);
        // printf("inv(A)*D = [%f %f %f]\n", dx, dy, dz);
        printf("X = [%f, %f, %f]\n", X[0], X[1], X[2]);
        assert(0);
    }
#endif

    if(Y!=NULL)
    { Y[0] = X[0]; Y[1] = X[1]; Y[2] = X[2]; }

    return dist;
}

/* Don't compile this when building object files */
#ifdef standalone

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
        double de = elli_pdist(E, p);
        double dg = elli_gdist(E, p, NULL);

        if( fabs(fabs(de) - fabs(r-norm3(p))) > 1e-6)
        {
            printf("elli_pdist = %f, |1-||p||| = %f\n", de, fabs(r-norm3(p)));
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
            elli_show(E);
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
    /* Random points for ellipse that has different axes
     * elli_gdist has it's own debugging method as long as
     * NDEBUG isn't defined
     */

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
    printf("    Ed = %f %f %f = E - %f\n", Ei->a, Ei->b, Ei->c, delta);
    printf("    Testing that f(p) < d for p \\ in Ed\n");



    E->maxiter = 1000;
//  elli_show(E);
//  elli_show(Ei);
    for(size_t kk = 0; kk<20000; kk++)
    {
        double X[3];
        double Y[3];
        X[0] = (double) rand()/ (double) RAND_MAX;
        X[1] = (double) rand()/ (double) RAND_MAX;
        X[2] = (double) rand()/ (double) RAND_MAX;
        elli_scale(Ei, X); // scale X to be on Ei
        double d = elli_gdistR(E, X, Y);
        if(d>=delta)
        {
            printf("No, for point %zu:", kk);
            printf("delta: %f d(xi, E) = %f\n", delta, eudist3(X,Y));
            assert(0);
        }
    }
    printf("    @");
    free(E);
    free(Ei);

}

int main(int argc, char ** argv)
{
#ifdef NDEBUG
    printf("Please turn on debugging or some checks will not be performed!\n");
#endif


    smalltest();

    if(argc == 4)
    {
        /* Test speed vs precision for a specific ellipse
         * N points on the ellipsoid are generated. Then they are
         * displaced along the normal and it is checked how well the
         * methods find back the original point.
         */
        size_t N = 1e7;
        double dist = 2*0.08; // distance from E
        size_t iterlast = 7; // For method with variable number of iterations
        double X[6];
        for(int kk = 1; kk<4; kk++)
        {
            X[kk-1] = atof(argv[kk]);
        }
        elli * E = elli_new(X[0], X[1], X[2]);

        E->eps = 1e-5;
        printf("eps: %e\n", E->eps);

        printf("Generating %zu random points on a sphere\n", N);
        double * P = getPointsOnSphere(N);
        assert(P != NULL);
        printf("Moving points from sphere to ellipse\n");
        for(size_t kk = 0; kk<N; kk++){
            elli_scale(E, P+3*kk);}

        printf("Generating reference points around the ellipse surface\n");
        double * Q = malloc(3*N*sizeof(double));
        assert(Q != NULL);
        for(size_t kk = 0; kk<N; kk++)
        {
            double normal[3];
            elli_nnormal(E, P+3*kk, normal);
            for(size_t ll = 0; ll<3; ll++)
            {
                if(1) // ll>N/2)
                {
                    Q[3*kk+ll] = P[3*kk+ll] - dist*normal[ll];
                } else {
                    Q[3*kk+ll] = P[3*kk+ll] + dist*normal[ll];
                }
            }
        }

        // Allocate space for output
        double * Y = malloc(3*N*sizeof(double));
        assert(Y != NULL);

        /* Now we have pairs (P_i, Q_i) such that P_i is the closest
         * point to Q_i on E
         * */
        struct timespec tstart;
        struct timespec tend;

        printf("Method    \tmaxerr    \trelerr    \tangerr    \tt       \ttpp\n");

        // 1X scaling
        clock_gettime(CLOCK_REALTIME, &tstart);
        double maxerror = 0;
        double maxnerror = 0;

        for(size_t kk = 0; kk<N; kk++)
        {
            memcpy(Y+3*kk, Q+3*kk, 3*sizeof(double));
            elli_scale(E, Y+3*kk);
        }
        clock_gettime(CLOCK_REALTIME, &tend);

        for(size_t kk = 0; kk<N; kk++)
        {
            double error = eudist3(P+3*kk, Y+3*kk);
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
        printf("scaling1");
        printf("\t%e\t%f", maxerror,100.0*maxerror/dist);
        printf("\t%e", maxnerror);
        double ttotal = clockdiff(&tstart, &tend);
        printf("\t%e\t%e\n", ttotal, N/ttotal);

        // Scaling x2
        clock_gettime(CLOCK_REALTIME, &tstart);
        maxerror = 0;
        maxnerror = 0;
        for(size_t kk = 0; kk<N; kk++)
        {
            memcpy(Y+3*kk, Q+3*kk, 3*sizeof(double));
            elli_scale2(E, Y+3*kk);
        }
        clock_gettime(CLOCK_REALTIME, &tend);

        for(size_t kk = 0; kk<N; kk++)
        {
            double error = eudist3(P+3*kk, Y+3*kk);
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
        printf("scaling2");
        printf("\t%e\t%e", maxerror,100.0*maxerror/dist);
        printf("\t%e", maxnerror);
        ttotal = clockdiff(&tstart, &tend);
        printf("\t%e\t%e\n", ttotal, N/ttotal);

        for(size_t ii = 1; ii<iterlast;  ii++)
        {
            if(ii + 1 == iterlast)
            { E->maxiter = 99; }
            else { E->maxiter = ii; }

            clock_gettime(CLOCK_REALTIME, &tstart);
            maxerror = 0;
            maxnerror = 0;
            for(size_t kk = 0; kk<N; kk++)
            {
                elli_gdist(E, Q+3*kk, Y+3*kk);
            }
            clock_gettime(CLOCK_REALTIME, &tend);
            for(size_t kk = 0; kk<N; kk++)
            {
                double error = eudist3(P+3*kk, Y+3*kk);
                if(error > maxerror)
                    maxerror = error;

                double nP[3];
                double nPE[3];
                elli_nnormal(E, P+3*kk, nP);
                elli_nnormal(E, Y+3*kk, nPE);
                double nerror = 360.0/M_PI*acos(sprod(nP, nPE));
                if(nerror > maxnerror)
                    maxnerror = nerror;

            }
            printf("bektas-%zu", E->maxiter);
            printf("\t%e\t%f", maxerror,100.0*maxerror/dist);
            printf("\t%e", maxnerror);
            ttotal = clockdiff(&tstart, &tend);
            printf("\t%e\t%e\n", ttotal, N/ttotal);
        }

        for(size_t ii = 1; ii<iterlast;  ii++)
        {
            if(ii + 1 == iterlast)
            { E->maxiter = 99; }
            else { E->maxiter = ii; }

            clock_gettime(CLOCK_REALTIME, &tstart);
            maxerror = 0;
            maxnerror = 0;
            for(size_t kk = 0; kk<N; kk++)
            {
                elli_gdistL(E, Q+3*kk, Y+3*kk);
            }
            clock_gettime(CLOCK_REALTIME, &tend);

            for(size_t kk = 0; kk<N; kk++)
            {
                double error = eudist3(P+3*kk, Y+3*kk);
                if(!isfinite(error))
                    error = 99e99;

                if(0)
                    printf("P=%f %f %f Y = %f %f %f\n",
                           P[3*kk], P[3*kk+1], P[3*kk+2],
                           Y[3*kk], Y[3*kk+1], Y[3*kk+2]);

                if(error > maxerror)
                    maxerror = error;

                double nP[3];
                double nPE[3];
                elli_nnormal(E, P+3*kk, nP);
                elli_nnormal(E, Y+3*kk, nPE);
                double nerror = 360.0/M_PI*acos(sprod(nP, nPE));
                if(nerror > maxnerror)
                    maxnerror = nerror;

            }
            printf("lagrange%zu", E->maxiter);
            printf("\t%e\t%f", maxerror,100.0*maxerror/dist);
            printf("\t%e", maxnerror);
            ttotal = clockdiff(&tstart, &tend);
            printf("\t%e\t%e\n", ttotal, N/ttotal);
        }

        for(size_t ii = 1; ii<iterlast;  ii++)
        {
            if(ii + 1 == iterlast)
            { E->maxiter = 99; }
            else { E->maxiter = ii; }


            clock_gettime(CLOCK_REALTIME, &tstart);
            maxerror = 0;
            maxnerror = 0;
            for(size_t kk = 0; kk<N; kk++)
            {
                elli_gdistR(E, Q+3*kk, Y+3*kk);
            }
            clock_gettime(CLOCK_REALTIME, &tend);

            for(size_t kk = 0; kk<N; kk++)
            {
                double error = eudist3(P+3*kk, Y+3*kk);
                if(!isfinite(error))
                    error = 99e99;

                if(0)
                    printf("P=%f %f %f Y = %f %f %f\n",
                           P[3*kk], P[3*kk+1], P[3*kk+2],
                           Y[3*kk], Y[3*kk+1], Y[3*kk+2]);

                if(error > maxerror)
                    maxerror = error;

                double nP[3];
                double nPE[3];
                elli_nnormal(E, P+3*kk, nP);
                elli_nnormal(E, Y+3*kk, nPE);
                double nerror = 360.0/M_PI*acos(sprod(nP, nPE));
                if(nerror > maxnerror)
                    maxnerror = nerror;
            }
            printf("polynomial%zu", E->maxiter);
            printf("\t%e\t%f", maxerror,100.0*maxerror/dist);
            printf("\t%e", maxnerror);
            ttotal = clockdiff(&tstart, &tend);
            printf("\t%e\t%e\n", ttotal, N/ttotal);
        }


        free(Q);
        free(E);
        free(Y);
        return 0;
    }


    if(argc == 7)
    {
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
        exit(1);
    }

    for(double r = .8; r<1.2; r=r+0.01)
    {
        compare_sphere(r);
    }

    geodesic_random();

    printf("All tests passed\n");
    return 0;
}

#endif
