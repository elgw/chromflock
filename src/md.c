#include "fast_prng/normal.h"
#define INLINED inline __attribute__((always_inline))

/* Molecular dynamics
*/

/* TODO:
 * - sort out the mess about the current steps,
 *   the steps should be enumerated both locally, at the current temperature
 *   as well as globally, over all temperatures. This is because some
 *   parameters depends on the global and some on the local step
 */

INLINED double eudist3(double * A, double * B)
{
  /* Euclidean distance between two 3D-vectors */
  return sqrt( pow(A[0]-B[0], 2) + pow(A[1]-B[1], 2) + pow(A[2]-B[2], 2));
}


INLINED double chrpd(double * X, optparam * p)
{
  // Monitor sum of pairwize distances of chr 1,
  // assuming it is correlated to compaction and volume
  // it should give an indication of structure clustering

  if(p->L != NULL)
  {

    size_t last = 0;
    for(size_t pp = 0; pp<p->N; pp++)
    {
      if(p->L[pp] == 1)
      {
        last = pp;
      }
    }

    double pd = 0;
    for(size_t kk = 0; kk<last; kk++)
    {
      for(size_t ll = kk+1; ll<last; ll++)
      {
        pd += eudist3(X+3*kk, X+3*ll);
      }
    }

    return pd;
  }
  return -1;
}

INLINED static void rand3d(double * d)
{
  // Provide a random 3d unit length vector
  double n = 2;
  do {
    d[0] = 2*(double) rand() / (double) RAND_MAX - 1;
    d[1] = 2*(double) rand() / (double) RAND_MAX - 1;
    d[2] = 2*(double) rand() / (double) RAND_MAX - 1;
    n = norm3(d);
  } while(n>1);
  d[0]/=n;
  d[1]/=n;
  d[2]/=n;
  return;
}

void comforce(optparam * restrict p, double * restrict X, double * restrict G)
{
  if(p->L == NULL)
  {
    printf("No L, can't compute com-force!\n");
    return;
  }

  // mX mean X positions
  size_t size_mX = 3*256*sizeof(double);
  double * mX = malloc(size_mX);
  memset(mX, 0, size_mX);
  // Number of points per label
  size_t * nL = malloc(256*sizeof(size_t));
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
  // Normalize by the number of points

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

  for(size_t kk = 0; kk < p->N; kk++)
  {
    uint8_t label = p->L[kk];
    double dist = eudist3(X+3*kk, mX+3*label);
    if(dist > 0)
    {
      for(size_t idx = 0; idx<3; idx++)
      {
        G[3*kk+idx] += 1*(X[3*kk+idx]-mX[3*label + idx])*pow(dist, 2) ;
      }
    }
  }

  free(nL);
  free(mX);

}


int dynamic(double * X, 
    optparam * p, 
    double Fb) // Brownian force
{

printf("Fb = %f\n", Fb);
  // TODO: 
  // - monitor errors and compaction to select good T curve. Use same seed or average over many
  //   structures.
  //   - Possibly different curve when the beads are initialized for the first time
  // - cap velocities?

  normal_setup();

  // Prepare settings
  size_t maxiter = p->maxiter;
  size_t nIterations = 0; // count

  size_t steps = 10000;
  if(steps>maxiter)
  {
    steps = maxiter;
  }

  conf fconf;
  fconf.r0 = p->r0;
  fconf.kVol = p->kVol;
  fconf.kSph = p->kSph;
  fconf.kInt = p->kInt; // Function of current iteration, see below
  fconf.kRad = p->kRad;
  fconf.dInteraction = 2.1*fconf.r0; // Function of current iteration, see below

  fconf.nIPairs = p->NI; // Only for err2

  // Initialization
  struct timespec ta, tb;
  clock_gettime(CLOCK_MONOTONIC, &ta);

  // Gradient
  double * g = malloc(p->N*3*sizeof(double)); // for gradient
  memset(g, 0, p->N*3*sizeof(double)); // zero initial velocity
  // Velocity
  double * v = malloc(p->N*3*sizeof(double)); // for velocity
  memset(v, 0, p->N*3*sizeof(double)); // zero initial velocity

  double * Xm = malloc(3*p->N*sizeof(double));
  // Set last position to be the current position, i.e. initial velocity will be 0
  for(size_t pp = 0; pp< 3*p->N ; pp++)
  {
    Xm[pp] = X[pp]; 
  }

  // see https://github.com/elgw/yeastSimulator/blob/master/sim.c
  int brownian = 1;

  // TODO: 
  // - shape of cooling off?
  // - stopping criteria?
  // - a period with no Fb at the end?

  double errmax = p->errstop; // 0.01;
  double gmax = p->gradstop; // 0.01;
  double error = 2*errmax;
  double gnorm = 2*gmax;
  double dt = 0.15;
  double damp = 0.5; // dampening (.8) faster convergence than .4

  do{
    for(size_t ss = 0; ss<steps; ss++)
    {
      int showStep = 0;

      if(ss % 1500 == 0 || ss == steps-1)
      {
        showStep = 1;
      }

      // 1. Estimate velocities
      // note, per component
      for(size_t pp = 0 ; pp < 3*p->N ; pp++)
      {
        v[pp] = (X[pp] - Xm[pp]) / (2 * dt);
      }

      // 2. Calculate Forces

      nIterations++;
      p->iter_final++;

      /* Alber style to avoid local minima */
      if(1)
      {
        double beta = 5.0;
        // WARNING: Assumes that 7 temperature cycles are used
        double n = (double) p->iter_final/ ((double) p->maxiter*7); // [0, 1]
        assert(n<=1);
        if(n>1)
        {
          n = 1;
        }

        fconf.kInt = p->kInt* 0.5*(1.0+erf((n-.5)*beta));
        if(fconf.kInt<0.1)
        {
          fconf.kInt = 0.1;
        }
        /* Define the contact range between beads. When beads are further away than the contact range
         * they will be attracted.
         * On page 64 of supplementary materials of nbt2057 it is stated that it goes from 1.1 to 2 times (r0+r1)
         * However, in PGS/alab/modeling.py the contact range is set to 2*(r0+r1) independent on iteration.
         * n \in [0,1], proportion of steps taken */
        fconf.dInteraction = (1+0.1+0.9*n)*2*fconf.r0;
        //printf("fconf.dInteraction: %f (%f r0)\n", fconf.dInteraction , fconf.dInteraction/fconf.r0);
      }

      /* Stevens style to avoid local minima */
      if(0) 
      {
        double beta = 5.0;
        // WARNING: Assumes that 7 temperature cycles are used
        double n = (double) p->iter_final/ ((double) p->maxiter*7);
        assert(n<=1);
        assert(n>0);      
        fconf.kVol = p->kVol * 0.5*(1.0+erf((n-.5)*beta));
      }

      // 2.1 Gradient of functional
      grad3(X, 
          p->N,
          p->R,
          p->I,
          g,
          &fconf);

      // Calculate gradient 2-norm
      if(showStep == 1)
      {
        gnorm = 0;
        for(size_t kk = 0; kk < 3*p->N; kk++)
        {
          gnorm += pow(v[kk], 2);
        }
        gnorm = sqrt(gnorm);
      }

      if(p->compress)
      {
        comforce(p, X, g);
      }

      // 2.2 Brownian force
      if(brownian)
      {
        // double Fb = (1-exp(-(double) steps/ ((double) ss+1)));
        // double Fb = 2*(1.0-((double) ss/ (double) steps)); // unstable
        //double Fb = (1.0-(((double) ss +0.1*steps)/ (double) steps));
        //Fb = .8; // !!! XXX 
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
            double Frnd = .10;
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
      }

      // 2.3 Dampening
      for(size_t kk = 0; kk < 3*p->N; kk++)
      {
        g[kk] = g[kk] + damp*v[kk]; 
      }

      if(0){ // Average gradient norm per bead
        double gavg = 0;
        for(size_t kk = 0; kk < p->N; kk++)
        {
          gavg += norm3(g+3*kk);
        }
        printf("Average gradient: %f\n", gavg/p->N);
      }

      // 3. Update X
      for(size_t pp = 0 ; pp < 3*p->N ; pp++)
      {
        double xt = X[pp];
        X[pp] = 2*X[pp] - Xm[pp] - g[pp]*pow(dt,2);
        Xm[pp] = xt; // Update Xm to reflect the previous X-value
      }

      if(showStep == 1)
      {
        double pd1 = chrpd(X, p);

        //radpos(X, p, stdout);

        error = err3(X,
            p->N,
            p->R,
            p->I,
            &fconf);
        printf("%zu, E: %f, ||G||: %f pd1: %f\n", ss, error, gnorm, pd1);

        // Abort if the error is small
        if(error < errmax)
        {
          printf("STOP!\n");
          ss = steps;
        }
        if(gnorm < gmax)
        {
          printf("STOP!\n");
          ss = steps;
        }

      }

      if(run == 0)
      { 
        // If Ctrl+C/interrupt was captured, abort.
        ss = steps;
      }

    }
  }  
  while( (nIterations < maxiter) 
      && (gnorm >=gmax)
      && (error >= errmax)
      && (run == 1));

  /* At final step, report back */
  double errorFinal = err3(X,
      p->N,
      p->R,
      p->I,
      &fconf);

  double gradFinal = 0;
  for(size_t kk = 0; kk < p->N; kk++)
  {
    gradFinal += pow(v[kk], 2);
  }
  gradFinal = sqrt(gradFinal);

  clock_gettime(CLOCK_MONOTONIC, &ta);

  p->err_final = errorFinal;
  p->grad_final = gradFinal;
  p->time_final = clockdiff(&ta, &tb);

  free(v);
  free(g);
  free(Xm);

  return 0;
}


int dynamic_old(double * X, 
    optparam * p, 
    double Fb) // Brownian force
{

  // TODO: 
  // - monitor errors and compaction to select good T curve. Use same seed or average over many
  //   structures.
  //   - Possibly different curve when the beads are initialized for the first time
  // - cap velocities?


  // Prepare settings
  size_t maxiter = p->maxiter;
  size_t nIterations = 0; // count

  size_t steps = 10000;
  if(steps>maxiter)
  {
    steps = maxiter;
  }

  conf fconf;
  fconf.r0 = p->r0;
  fconf.kVol = p->kVol;
  fconf.kSph = p->kSph;
  fconf.kInt = p->kInt;
  fconf.kRad = p->kRad;

  fconf.dInteraction = 2.1*fconf.r0;
  fconf.nIPairs = p->NI; // Only for err2

  // Initialization
  struct timespec ta, tb;
  clock_gettime(CLOCK_MONOTONIC, &ta);

  // Gradient
  double * g = malloc(p->N*3*sizeof(double)); // for gradient
  memset(g, 0, p->N*3*sizeof(double)); // zero initial velocity
  // Velocity
  double * v = malloc(p->N*3*sizeof(double)); // for velocity
  memset(v, 0, p->N*3*sizeof(double)); // zero initial velocity

  // see https://github.com/elgw/yeastSimulator/blob/master/sim.c
  int brownian = 1;

  // TODO: 
  // - shape of cooling off?
  // - stopping criteria?
  // - a period with no Fb at the end?

  double errmax = p->errstop; // 0.01;
  double gmax = p->gradstop; // 0.01;
  double error = 2*errmax;
  double gnorm = 2*gmax;
  do{
    for(size_t ss = 0; ss<steps; ss++)
    {
      int showStep = 0;

      if(ss % 500 == 0 || ss == steps-1)
      {
        showStep = 1;
      }

      nIterations++;
      p->iter_final++;

      double kInt = p->kInt;
      p->kInt = (1-Fb)*(1-Fb)*kInt;

      grad3(X, 
          p->N,
          p->R,
          p->I,
          g,
          &fconf);

      // Calculate gradient 2-norm
      if(showStep == 1)
      {
        gnorm = 0;
        for(size_t kk = 0; kk < 3*p->N; kk++)
        {
          gnorm += pow(v[kk], 2);
        }
        gnorm = sqrt(gnorm);
      }

      if(p->compress)
      {
        comforce(p, X, g);
      }


      if(brownian)
      {
        // double Fb = (1-exp(-(double) steps/ ((double) ss+1)));
        // double Fb = 2*(1.0-((double) ss/ (double) steps)); // unstable
        //double Fb = (1.0-(((double) ss +0.1*steps)/ (double) steps));
        //Fb = .8; // !!! XXX 
        if(Fb<0)
          Fb = 0;

        if(Fb>0)
        {
          for(size_t kk = 0; kk < p->N; kk++)
          {
            double d[] = {0,0,0};
            rand3d(&d[0]);
            double Frnd = 0.5*((double) rand() / (double) RAND_MAX - 0.5);
            for(size_t idx =0 ; idx<3; idx++)
            {
              g[3*kk+idx] += Fb*Frnd*d[idx];
              // g[3*kk+idx] = (1-Fb)*g[3*kk+idx] + Fb*Frnd*d[idx];

            }
          }
        }
      }

      p->kInt = kInt;

      //      double damp = 0.8; // dampening (.8) faster convergence than .4
      double damp = 0.8; // dampening (.8) faster convergence than .4
      double deltat = 5; // 1=ok
      for(size_t kk = 0; kk < 3*p->N; kk++)
      {
        // e-3 -> 0.3 e-2 -> 0.03 e-1 -> 0.01
        v[kk] = damp*v[kk] - 1e-2*g[kk]; 
        X[kk] = X[kk] + deltat*v[kk];
      }


      if(showStep == 1)
      {
        double pd1 = chrpd(X, p);

        //radpos(X, p, stdout);

        error = err3(X,
            p->N,
            p->R,
            p->I,
            &fconf);
        printf("%zu, E: %f, ||G||: %f pd1: %f\n", ss, error, gnorm, pd1);

        // Abort if the error is small
        if(error < errmax)
        {
          printf("STOP!\n");
          ss = steps;
        }
        if(gnorm < gmax)
        {
          printf("STOP!\n");
          ss = steps;
        }

      }
      if(run == 0)
      {
        // Abort!
        ss = steps;
      }

    }
  }  
  while( (nIterations < maxiter) 
      && (gnorm >=gmax)
      && (error >= errmax)
      && (run == 1));

  /* At final step, report back */
  double errorFinal = err3(X,
      p->N,
      p->R,
      p->I,
      &fconf);

  double gradFinal = 0;
  for(size_t kk = 0; kk < p->N; kk++)
  {
    gradFinal += pow(v[kk], 2);
  }
  gradFinal = sqrt(gradFinal);

  clock_gettime(CLOCK_MONOTONIC, &ta);

  p->err_final = errorFinal;
  p->grad_final = gradFinal;
  p->time_final = clockdiff(&ta, &tb);

  free(v);
  free(g);

  return 0;
}

