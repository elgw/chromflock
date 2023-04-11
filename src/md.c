/* Molecular dynamics
*/

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

INLINED double chrpd(const double * X, const optparam * p)
{
    /* Monitor sum of pairwize distances of chr 1,
     *  assuming it is correlated to compaction and volume
     * it should give an indication of structure clustering */

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
    /* Provide a random 3d unit length vector */
  double n = 2;
  do {
    d[0] = 2*(double) rand() / (double) RAND_MAX - 1;
    d[1] = 2*(double) rand() / (double) RAND_MAX - 1;
    d[2] = 2*(double) rand() / (double) RAND_MAX - 1;
    n = norm3(d);
  } while( n>1 || n < 1e-6);
  /* If the point is in the unit sphere,
     and not too close to origo, continue */
  d[0]/=n;
  d[1]/=n;
  d[2]/=n;
  return;
}

/* Apply a force that attracts each bead to the centre of
   mass of it's chromosome */
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


int dynamic(double * X,
    optparam * p,
            double Fb) /* Brownian force */
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
  double * g = malloc(p->N*3*sizeof(double));
  memset(g, 0, p->N*3*sizeof(double)); /* zero initial velocity */
  // Velocity
  double * v = malloc(p->N*3*sizeof(double));
  memset(v, 0, p->N*3*sizeof(double)); /* zero initial velocity */
  double * Xm = malloc(3*p->N*sizeof(double));
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
        /* Alber style heuristics */
      double beta = 5.0;
      /* WARNING: Assumes that 7 temperature cycles are used */
      double n = (double) p->iter_final/ ((double) p->maxiter*7); // [0, 1]
      assert(n<=1);
      if(n>1)
      {
        n = 1;
      }

      fconf.kInt = p->kInt* 0.5*(1.0+erf((n-.5)*beta));
      if(fconf.kInt<0.1*p->kInt)
      {
        fconf.kInt = 0.1*p->kInt;
      }
      /* Define the contact range between beads. When beads are further away than the contact range
       * they will be attracted.
       * On page 64 of supplementary materials of nbt2057 it is stated that it goes from 1.1 to 2 times (r0+r1)
       * However, in PGS/alab/modeling.py the contact range is set to 2*(r0+r1) independent on iteration.
       * n \in [0,1], proportion of steps taken */
      fconf.dInteraction = (1+0.1+0.9*n)*2*fconf.r0;
      //printf("fconf.dInteraction: %f (%f r0)\n", fconf.dInteraction , fconf.dInteraction/fconf.r0);

    } else { /* Settings from lua script */
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
      v[pp] = (X[pp] - Xm[pp]) / (2 * dt);
    }

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

    /* 3. Update X */
    for(size_t pp = 0 ; pp < 3*p->N ; pp++)
    {
      double xt = X[pp];
      X[pp] = 2*X[pp] - Xm[pp] - g[pp]*pow(dt,2);
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

      //      double pd1 = chrpd(X, p);

      error = err3(X,
          p->N,
          p->R,
          p->I,
          &fconf);
      logwrite(p, 2, "    Iter: %6zu, E: %e, ||G||: %e\n", iter, error, gnorm);
    }

  } while( (iter < maxiter) && (run == 1) && (luaquit == 0));

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
