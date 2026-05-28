#include "ellipsoid.h"

static double
eudist3(const double * a, const double * b)
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

static double
norm3(const double * a)
{
    return sqrt(
        pow(a[0], 2) +
        pow(a[1], 2) +
        pow(a[2], 2));
}

elli * elli_new(double a, double b, double c)
{
    elli * E = calloc(1, sizeof(elli));
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

void elli_print(FILE * fid, elli * E)
{
    fprintf(fid, "Ellipse: x^2/a^2 + y^2/b^2 + z^2/c^2 = 1, maxiter: %zu, eps: %e\n",
           E->maxiter, E->eps);
    fprintf(fid, "         a=%f, b=%f, c=%f, volume=%f\n",
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

double
elli_radius(const elli * E, const double * p)
{

    double t = sqrt( 1.0/(  pow(p[0]/E->a, 2)
                            + pow(p[1]/E->b, 2)
                            + pow(p[2]/E->c, 2)
                         ));

    return 1.0/t;
}

double
elli_dist_p(const elli * E,
           const double * restrict p,
           double * restrict X)
{

    double _ip[] = {0,0,0};
    double * ip;
    if(X == NULL) {
        ip = _ip;
    } else {
        ip = X;
    }
    memcpy(ip, p, 3*sizeof(double));
    elli_project(E, ip);

    return eudist3(p, ip);
}

double
elli_getScale(const elli * restrict E, const double * restrict X)
{
    return sqrt(elli_getScale2(E, X));
}

double
elli_getScale2(const elli * restrict E, const double * restrict X)
{
    return pow(X[0]/E->a,2) + pow(X[1]/E->b,2) + pow(X[2]/E->c,2);
}

void
elli_project(const elli * E, double * X)
{
    double z = sqrt(pow(X[0]/E->a,2) + pow(X[1]/E->b,2) + pow(X[2]/E->c,2));
    if(z == 0)
    {
        X[0] = 0; X[1] = 0; X[2] = 0;
    } else {
        X[0]*=1.0/z; X[1]*=1.0/z; X[2]*=1.0/z;
    }
}


void elli_normal(const elli * E,
                 const double * restrict X,
                 double * restrict N)
{
    N[0] = 2.0*X[0]/pow(E->a,2);
    N[1] = 2.0*X[1]/pow(E->b,2);
    N[2] = 2.0*X[2]/pow(E->c,2);
}

void
elli_nnormal(const elli * E,
             const double * restrict X,
             double * restrict N)
{
    elli_normal(E, X, N);
    double n = norm3(N);
    for(int kk = 0; kk<3; kk++)
    {
        N[kk]/=n;
    }
}

static double
f_div_df(const double l,
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

static double
elli_dgistR_fun(const double * E, const double * P, double alpha)
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

double
elli_gdistR(const elli * ellipse,
            const double * restrict P,
            double * restrict X)
{

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

double
elli_gdistL(const elli * ellipse,
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
    if(norm3(P) > norm3(Y))
    {
        return eudist3(P, Y);
    } else {
        return -eudist3(P,Y);
    }
}

double
elli_gdist(const elli * ellipse,
           const double * restrict P,
           double * restrict Y)
{
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

    if(Y != NULL) {
        Y[0] = X[0];
        Y[1] = X[1];
        Y[2] = X[2];
    }

    return dist;
}

double
elli_gdistS(const elli * E,
            const double * P)
{
    double x2 = pow(P[0], 2);
    double y2 = pow(P[1], 2);
    double z2 = pow(P[2], 2);

    // F(x)
    double num = x2/E->a2 + y2/E->b2 + z2/E->c2 - 1.0;
    // ||grad F(x)||
    double den = 2*sqrt(x2/(E->a2*E->a2)
                        + y2/(E->b2*E->b2)
                        + z2/(E->c2*E->c2));
    return - num/den;
}
