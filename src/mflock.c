#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <unistd.h>
#include <signal.h>
#include <stdarg.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <math.h>

#include "mflock.h"
#include "cmmwrite.h"
#include "functional.h"
#include "wio.h"

#ifdef SDL
#include <pthread.h>
#include "liveview.h"
#include "liveview.c"
#endif

/* Ctrl+C switches run to 0 which causes
 * the solver to stop asap */

#define INLINED inline __attribute__((always_inline))

void errline(int line)
{
  printf("L%d\n", line);
  fflush(stdout);
  return;
}

static volatile int run = 1;
static void stoprun(int ignore)
{
  run = 0;
}


static double clockdiff(struct timespec* start, struct timespec * finish)
{

  double elapsed = (finish->tv_sec - start->tv_sec);
  elapsed += (finish->tv_nsec - start->tv_nsec) / 1000000000.0;
  return elapsed;
}


static double norm3(double * X)
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

void radpos(double * X, optparam * p, FILE * f)
{
  // Show radial positions of each chromosome, using centre of mass
  // an alternative would be to show mean radius.
  //
  if(p->L == NULL)
  {
    printf("No L!\n");
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


#include "md.c" // Molecular dynamics, provides: dynamic


void param_summary(optparam * p, double * X, FILE * f)
{
  assert(p->N>0);

  fprintf(f, "\n");
  fprintf(f, ">> Optimization summary:\n");
  fprintf(f, "   Final iterations: %zu\n", p->iter_final);
  fprintf(f, "   Total time: %zu s\n", p->time_final);
  fprintf(f, "   Final gradient: %e\n", p->grad_final);
  fprintf(f, "   Final error: %e\n", p->err_final);

  // X: mean, max, min
  double mex = 0, mey = 0, mez = 0;
  double md = 10e99;
  double mix = md, miy = md, miz = md;
  double max = -md, may = -md, maz = -md;
  double mer = 0, mir = 10e99, mar = 0;

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

#ifndef DNDEBUG
  fprintf(f, ">> Structure summary:\n");
  fprintf(f, "              X       Y       Z       R\n");
  fprintf(f, "   Max:  % .3f, % .3f, % .3f, % .3f\n", max, may, maz, mar);
  fprintf(f, "   Mean: % .3f, % .3f, % .3f, % .3f\n", mex, mey, mez, mer);
  fprintf(f, "   Min:  % .3f, % .3f, % .3f, % .3f\n", mix, miy, miz, mir);
#endif

  if(run == 0)
  {
    fprintf(f, "WARNING: abnormal exit (Ctrl+c was pressed?)\n");
  }

}

void param_readW(optparam * p)
{

  loginfo(p, "Reading pairwise interactions from %s\n", p->wfname);

  size_t fsize = 0;
  p->W = wio_read(p->wfname, &fsize);
  if(fsize == 0)
  {
    printf("Failed reading file!\n");
    exit(1);
  }

  loginfo(p, " W = [ %d, %d, %d, %d, ..., %d\n", 
      (int) p->W[0], 
      (int) p->W[1], 
      (int) p->W[2], 
      (int) p->W[3], 
      p->W[fsize-1]);

  // Set the number of elements
  p->N = (size_t) sqrt(fsize);
  assert(p->N*p->N == fsize);
  loginfo(p, " %zu points\n", p->N);

  // Initiate the list with pairs in contact
  // I.e. p->I and p->NI

  // Count them
  size_t nPairs = 0;
  for(size_t kk = 0; kk<p->N; kk++)
  {
    for(size_t ll = kk+1; ll<p->N; ll++)
    {
      if(p->W[kk+ll*p->N] == 1)
      {
        nPairs++;
      }
    }
  }

  printf(" Found %zu pairs\n", nPairs);

  if(nPairs*20 > p->N)
  {
    printf("WARNING: more than 20 pairs per bead!\n");
  }

  p->NI = nPairs;
  p->I = malloc(nPairs*2*sizeof(uint32_t));
  assert(p->I != NULL);

  // Construct the list
  size_t writepos = 0;
  for(size_t kk = 0; kk<p->N; kk++)
    for(size_t ll = kk+1; ll<p->N; ll++)
    {
      if(p->W[kk+ll*p->N] == 1)
      {
        p->I[writepos++] = kk;
        p->I[writepos++] = ll;
      }
    }

  return;
}

int param_readL(optparam * p)
{
  /* Read label matrix pointed to by p->lfname
  */

  loginfo(p, "Reading L-labels from %s\n", p->lfname);
  int readAsText = 0;
  p->L = malloc(p->N*sizeof(double));

  // Try to read as binary
  FILE * f = fopen(p->lfname, "rb");
  fseek(f, 0, SEEK_END); // seek to end of file
  size_t fsize = ftell(f); // get current file pointer
  fseek(f, 0, SEEK_SET); // seek back to beginning of file

  loginfo(p, "As uint8_t, file contains %zu numbers (%zu bytes)\n",
      fsize/sizeof(uint8_t), fsize);

  if(fsize/sizeof(uint8_t) == p->N)
  {
    size_t status = fread(p->L, sizeof(uint8_t), p->N, f);    
    fclose(f);
    if(status != p->N)
    {
      printf("Problems at %d\n", __LINE__);
    }
  }
  else
  {
    fclose(f);
    printf("! Wrong number if numbers in L!\n");
    printf("! Trying to read as text\n");
    readAsText = 1;    
  }

  if(readAsText == 1)
  {
    loginfo(p, "Reading %s as text, one float per line\n", p->lfname);
    f = fopen(p->lfname, "r");
    size_t nbuf = 1024;
    char * buf = malloc(nbuf*sizeof(char));

    size_t pos = 0;
    while (fgets(buf, nbuf, f) != NULL && (pos < p->N))
    {  
      p->L[pos++] = atof(buf);
    }
    fclose(f);
    free(buf);
  }

  logdebug(p, "L = [%u, %u, ..., %u]\n", p->L[0], p->L[1], p->L[p->N-1]);

  return 0;
}

int param_readR(optparam * p)
{
  /* Read GPSeq radius values as binary double.
   * If that does not work, try as text, one value per line
   */
  loginfo(p, "Reading R-values from %s\n", p->rfname);
  size_t nbytes = 0;
  p->R = (double *) wio_read(p->rfname, &nbytes);
  printf("Read %zu doubles\n", nbytes/sizeof(double));
  if(nbytes/sizeof(double) != p->N)
  {
    printf("Seems like %s is corrupt or does not exist\n", p->rfname);
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

  /*
  int readAsText = 0;
  p->R = malloc(p->N*sizeof(double));

  // Try to read as binary
  FILE * f = fopen(p->rfname, "rb");
  assert(f != NULL);

  fseek(f, 0, SEEK_END); // seek to end of file
  size_t fsize = ftell(f); // get current file pointer
  fseek(f, 0, SEEK_SET); // seek back to beginning of file

  loginfo(p, "As 64 bit double, file contains %zu numbers (%zu bytes)\n",
      fsize/8, fsize);

  if(fsize/8 == p->N)
  {
    size_t status = fread(p->R, 8, p->N, f);    
    fclose(f);
    if(status != p->N)
    {
      printf("Problems at %d\n", __LINE__);
    }
  }
  else
  {
    fclose(f);
    printf("! Wrong number if numbers in R!\n");
    printf("! Trying to read as text\n");
    readAsText = 1;    
  }

  if(readAsText == 1)
  {
    loginfo(p, "Reading %s as text, one float per line\n", p->rfname);
    f = fopen(p->rfname, "r");
    size_t nbuf = 1024;
    char * buf = malloc(nbuf*sizeof(char));

    size_t pos = 0;
    while (fgets(buf, nbuf, f) != NULL && (pos < p->N))
    {  
      p->R[pos++] = atof(buf);
    }
    fclose(f);
    free(buf);
  }

  logdebug(p, "R = [%f, %f, ..., %f]\n", p->R[0], p->R[1], p->R[p->N-1]);

  size_t nInf = 0;
  for(size_t kk = 0; kk<p->N; kk++)
  {
    if(!isfinite(p->R[kk]))        
    {
      nInf++;
    }
  }

  printf("%s contains %zu non-finite values which will be ignored.\n", p->rfname, nInf);

*/
  return 0;
}


void param_dumpX(optparam * p, double * X)
{
  size_t N = p->N;


  int write_R = 0;
  if(p->R != NULL)
    write_R = 1;

  if(run == 1)
  {    
    loginfo(p, "Writing final structure to %s\n", p->xoutfname);
  }
  else {
    p->xoutfname = realloc(p->xoutfname, sizeof(char)*(strlen(p->xoutfname)+5));
    strcat(p->xoutfname, ".0");
    //    sprintf(p->xoutfname, "%s.0", p->xoutfname);

    loginfo(p, "Writing non-finished structure to: %s\n", p->xoutfname);
  }

  loginfo(p, "Columns: x, y, z, r");
  if(write_R)
    loginfo(p, ", R");
  loginfo(p, "\n");

  FILE * f = fopen(p->xoutfname, "w");

  for(size_t kk = 0; kk<N; kk++)
  {
    fprintf(f, "%f, %f, %f, %f", X[3*kk], X[3*kk+1], X[3*kk+2], norm3(X+3*kk));    
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

  fprintf(f, "\n");
  fprintf(f, "### Parameters:\n");
  fprintf(f, "# Problem size: %zu points (%zu variables)\n", p->N, 3*p->N);
  fprintf(f, "## Model parameters\n");
  fprintf(f, "# Spring constant volume exclusion (repulsion): %f\n", p->kVol);
  fprintf(f, "# Spring constant volume interactions: %f\n", p->kInt);
  fprintf(f, "# Spring constant sphere confinement: %f\n", p->kSph);
  fprintf(f, "# Spring constant radial positioning: %f\n", p->kRad);
  fprintf(f, "# Bead radius %f (vol. occ. %f)\n", p->r0, volocc);
  fprintf(f, "## Optimization parameters\n");
  if(p->dynamic == 1)
  {
    fprintf(f, "# Optimization: simulated annealing / molecular dynamics\n");
    if(p->compress == 1)
    {
      fprintf(f, "# Compression: on\n");
    } else {
      fprintf(f, "# Compression: off\n");
    }
  } else {
    fprintf(f, "# Optimization: bfgs\n");
  }
  fprintf(f, "# Error stop: %e\n", p->errstop);
  fprintf(f, "# Gradient stop: %e\n", p->gradstop);
  fprintf(f, "# Line search tol: %e\n", p->linestop);
  fprintf(f, "# Max iterations: %zu\n", p->maxiter);
  fprintf(f, "# Max time: %zu (s)\n", p->maxtime);
  fprintf(f, "# random seed: %zu\n", p->rseed);
  fprintf(f, "# write compresseed cmm: %d\n", p->cmmz);

  if(p->wfname == NULL)
  {
    fprintf(f, "# W file: -not specified-  REQUIRED\n");
  } else
  {
    fprintf(f, "# W file: %s\n", p->wfname);
  }

  if(p->xfname == NULL)
  {
    fprintf(f, "# X file: -not specified-\n");
  } else
  {
    fprintf(f, "# X file: %s\n", p->xfname);
  }

  if(p->lfname == NULL)
  {
    fprintf(f, "# L file: -not specified-\n");
  } else
  {
    fprintf(f, "# L file: %s\n", p->lfname);
  }

  if(p->rfname == NULL)
  {
    fprintf(f, "# R file: -not specified-\n");
  } else
  {
    fprintf(f, "# R file: %s\n", p->rfname);
  }

  if(p->ofoldername == NULL)
  {
    fprintf(f, "# output folder not specified\n");
  }
  else
  {
    fprintf(f, "# output folder: %s\n", p->ofoldername);
  }

  fprintf(f, "# verbose level: %d\n", p->verbose);
  fprintf(f, "\n");
}

void usage()
{
  printf("Accepted arguments:\n");
  printf("Required:\n");
  printf(" -w <file>, --wFile <file>\n\tfile with contacts marked out (square)\n");
  printf(" -L <file>, --lFile <file>\b\tlabel matrix (uint8_t)\n");
  printf("Optional:\n");
  printf(" -x <file>, --xFile <file>\n\tfile with initial coordinates\n");
  printf(" -s N, --seed N\n\tseed for random number generator (defaults to random)\n");
  printf(" -r <file>, --rFile <file>\n\tfile with wanted radii (in conjunction with -G)\n");
  printf("\tall files should be encoded as 64-bit floats\n");
  printf(" -e errstop X, --errstop X\n\t stop condition on error\n");
  printf(" -g gstop, --gradstop gstop\n\tstop condition on gradient magnitude\n");
  printf(" -l ltol, --linestop ltol\n\ttolerance for the line search, 0=exact\n");
  printf(" -n N, --maxiter N\n\tmaximum number of iterations\n");
  printf(" -t N, --maxtme N\n\t time budget in seconds\n");

  // Model parameters
  printf(" -V kVol, --kVol kVol\n\tvolume exclusions spring constant\n");
  printf(" -I kInt, --kInt kInt\n\tinteraction spring constant\n");
  printf(" -S kSph, --kSph kSph\n\tsphere confinement spring constant\n");
  printf(" -G kRad, --kRad kRad\n\tgenome positioning spring constant\n");
  printf(" -R r0, --radius r0\n\tbead radius (sphere has radius 1)\n");

  printf(" -v level, --verbose level\n\tverbosity lowest=0\n");
  printf(" -o name, --outFolder name\n\twhere to store results\n");
  printf(" -z, --cmmz\n\twrite compressed cmm files (gzip)\n");
  printf(" -D, --dynamic\n\tuse molecular dynamics\n");
  printf(" -c, --compress\n\tenable chromosome compression (only for -D)\n");
  printf(" -a, --live\n\tenable live monitoring (only when compiled with SDL)\n");
  printf(" -h, --help\n\tshow this help message. For more info see 'man mflock'\n");
  printf(" -d, --defaults\n\tshow default settings for the parameters\n");


  printf("\n");
}

int argparsing(optparam * p, int argc, char ** argv)
{

  struct option longopts[] = {
    {"version",     no_argument,       NULL,   'i' }, 
    { "help",         no_argument,       NULL,   'h' },
    // Data
    { "wFile",        required_argument, NULL,   'w' },
    { "xFile",        required_argument, NULL,   'x' },
    { "rFile",        required_argument, NULL,   'r' },
    { "lFile",        required_argument, NULL,   'L' },
    { "outFolder",    required_argument, NULL,   'o' },
    // Settings
    { "gradstop",     required_argument, NULL,   'g' },
    { "errstop",      required_argument, NULL,   'e' },
    { "linestop",     required_argument, NULL,   'l' },
    { "maxiter",      required_argument, NULL,   'n' },
    { "maxtime",      required_argument, NULL,   't' },
    { "seed",         required_argument, NULL,   's' },
    { "verbose",      required_argument, NULL,   'v' },
    { "live",         no_argument, NULL,   'a' },
    { "cmmz",         no_argument, NULL,   'z' },
    // Settings / Forces
    { "kVol",        required_argument, NULL,   'V' },
    { "kSph",        required_argument, NULL,   'S' },
    { "kInt",        required_argument, NULL,   'I' },
    { "kRad",        required_argument, NULL,   'G' },
    // Settings / Geometry
    { "radius",        required_argument, NULL,   'R' },
    // Settings / Optimization
    { "dynamic",        no_argument,       NULL,   'D' },
    { "compress",       no_argument,       NULL,   'c' },
    { "defaults",       no_argument,       NULL,   'd' },
    { NULL,           0,                 NULL,   0   }
  };


  int ch;
  while((ch = getopt_long(argc, argv, "w:x:r:g:l:n:t:V:I:S:R:G:v:o:hMs:De:L:zcad", longopts, NULL)) != -1)
  {
    switch(ch) {
      case 'i':
        printf("Build date: %s, %s\n", __DATE__, __TIME__);
        printf("GIT HASH: %s\n", GIT_VERSION);
        printf("Compiler: %s\n", CC_VERSION);
        exit(0);
      case 'd':
        printf("Defaults:\n");
        param_show(p, stdout);
        exit(0);
        break;
      case 'w':
        p->wfname = malloc(strlen(optarg)+1);
        strcpy(p->wfname, optarg);
        break;
      case 'x':
        p->xfname = malloc(strlen(optarg)+1);
        strcpy(p->xfname, optarg);
        break;
      case 'r':
        p->rfname = malloc(strlen(optarg)+1);
        strcpy(p->rfname, optarg);
        break;
      case 'L':
        p->lfname = malloc(strlen(optarg)+1);
        strcpy(p->lfname, optarg);
        break;
      case 'g':
        p->gradstop = atof(optarg);
        break;
      case 'e':
        p->errstop = atof(optarg);
        break;
      case 'l':
        p->linestop = atof(optarg);
        break;
      case 'n':
        p->maxiter = atol(optarg);
        break;
      case 't':
        p->maxtime = atol(optarg);
        break;
      case 'V':
        p->kVol = atof(optarg);
        break;
      case 'S':
        p->kSph = atol(optarg);
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
      case 's':
        p->rseed = atol(optarg);
        break;
      case 'v':
        p->verbose = atoi(optarg);
        break;
      case 'o':
        p->ofoldername = malloc(strlen(optarg)+1);
        strcpy(p->ofoldername, optarg);
        break;
      case 'h':
        return(1);
      case 'z':
        p->cmmz = 1;
        break;
      case 'D':
        p->dynamic = 1;
        printf("Dynamic!\n");
        break;
      case 'c':
        p->compress = 1;
        break;
      case 'a':
        p->liveView = 1;
        break;
      default:
        return(1);        
    }
  }

  if(p->wfname == NULL)
  {
    printf("'W' file not set (-w) \n");
    printf("Generate one in MATLAB with:\n");
    printf("A = zeros(150,150)\n");
    printf("A(2:size(A,1)+1:end) = 1;\n");
    printf("A = A + A';\n");
    printf("A8 = uint8(A);\n");
    printf("fout = fopen('A150.dat', 'wb');\n");
    printf("fwrite(fout, A, 'uint8');\n");
    printf("fclose(fout)\n");
    printf("\n");

    return 1;
  }

  return 0;
}

optparam *  param_alloc(void)
{
  optparam * p = malloc(sizeof(optparam));

  struct timespec ts;
  //  timespec_get(&ts, TIME_UTC);
  clock_gettime(CLOCK_REALTIME, &ts);

  p->rseed = time(NULL)*getpid()*ts.tv_nsec;

  p->errstop = 0.05; // stop if error is less than this
  p->gradstop = 1e-5; // gradient norm
  p->linestop = .1; // precision in line search
  p->maxiter = 100000; // iterations
  p->maxtime = 60*60*10; // seconds

  p->grad_final = 0;
  p->iter_final = 0;
  p->time_final = 0;

  p->kVol = 1;
  p->kInt = 1;
  p->kSph = 1;
  p->kRad = 0;
  p->r0 = 0.03;

  p->N = 0;
  p->W = NULL;
  p->I = NULL;
  p->NI = 0;
  p->R = NULL;
  p->verbose = 1;
  p->dynamic = 0;
  p->compress = 0;
  p->cmmz = 0;
  p->liveView = 0;

  p->L = NULL;

  p->rfname = NULL;
  p->xfname = NULL;
  p->wfname = NULL;
  p->ofoldername = NULL;
  p->lfname = NULL;

  p->logf = NULL;

  // Create a suggestion for the output folder
  if(p->ofoldername == NULL)
  {

    p->ofoldername = malloc(1024*sizeof(char));

    size_t fn = 1;
    int okfolder = 0;

    struct stat st = {0};

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

int param_init(optparam * p)
{

  // Create output folder 
  loginfo(p, "Suggested output folder: %s\n", p->ofoldername);
  struct stat st = {0};
  if(stat(p->ofoldername, &st) == -1)
  {
    loginfo(p, "Creating output folder: %s\n", p->ofoldername);
    mkdir(p->ofoldername, 0770);
  } else {
    loginfo(p, "Warning: writing to existing folder\n");
  }


  // Set names of output files
  p->xoutfname = malloc(1024*sizeof(char));
  sprintf(p->xoutfname, "%s%s", p->ofoldername, "coords.csv");


  p->logfname = malloc(1024*sizeof(char));
  sprintf(p->logfname, "%s%s", p->ofoldername, "log.txt");


  // Open log
  p->logf = fopen(p->logfname, "w");
  assert(p->logf != NULL);

  if(p->logfname == NULL)
  {
    printf("Failed to open log file\n");
    return 1;
  }

  return 0;
}


void param_free(optparam * p)
{

  if(p->R != NULL)
    free(p->R);

  if(p->I != NULL)
    free(p->I);

  if(p->W != NULL)
    free(p->W);

  if(p->L != NULL)
    free(p->L);

  if(p->wfname != NULL)
    free(p->wfname);

  if(p->lfname != NULL)
    free(p->lfname);

  if(p->rfname != NULL)
    free(p->rfname);

  if(p->xfname != NULL)
    free(p->xfname);

  if(p->xoutfname != NULL)
    free(p->xoutfname);

  if(p->ofoldername != NULL)
    free(p->ofoldername);

  if(p->logfname != NULL)
    free(p->logfname);

  free(p);
  return;
}

void param_validate(optparam * p)
{
  if(p->W == NULL)
  {
    fprintf(stderr, "p->W == NULL\n");
    fprintf(stderr, "This indicates that a W file wasn't specified or that it could not be read\n");
    exit(1);
  }

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
}


void logwarn(optparam * p, const char *fmt, ...)
{
  if(p->verbose > 0)
  {
    va_list args, args2;
    va_start(args, fmt);
    va_copy(args2, args);
    // Prepend something like DBG:
    fprintf(stdout, "WARNING:");
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

void loginfo(optparam * p, const char *fmt, ...)
{
  if(p->verbose > 0)
  {
    va_list args, args2;
    va_start(args, fmt);
    va_copy(args2, args);
    // Prepend something like DBG:
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

void logdebug(optparam * p, const char *fmt, ...)
{
  if(p->verbose > 1)
  {
    va_list args, args2;
    va_start(args, fmt);
    va_copy(args2, args);
    // Prepend something like DBG:
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

int main(int argc, char ** argv)
{
  time_t starttime, nowtime;
  time(&starttime);

  optparam * p = param_alloc(); // Set default options

  if( argparsing(p, argc, argv)) 
  {
    usage();
    param_free(p); 
    return 1; 
  }

  param_init(p);

  param_readW(p);

  logdebug(p, "First value: %d, # values: %zu\n", (int) p->W[0], p->N);

  if(p->rfname != NULL)
  {
    param_readR(p);
  }

  if(p->lfname != NULL)
  {
    param_readL(p);
  }

  double * X = NULL;

  if(p->xfname != NULL)
  {
    X = malloc(3*p->N*sizeof(double));
    if(param_readX(p, X) != 0)
    {
      logdebug(p, "Could not open x-file, starting from random\n");
      free(X);
      X=NULL;
    }
  }

  if(X==NULL)
  {
    loginfo(p, "Using random initialization for X\n");
    srand(p->rseed);
    X = malloc(3*p->N*sizeof(double));
    for(size_t kk = 0; kk< p->N; kk++)
    {
      int accepted = 0;
      while(accepted == 0)
      {
        for(int idx =0; idx<3; idx++)
        {
          X[3*kk+idx] = 2*(rand()/(double) RAND_MAX-.5);    
        }
        if(norm3(X+3*kk)<1)
        { accepted = 1;}
      }
    }

    loginfo(p, "X[0] = %f\n", X[0]);
  }

  if(X == NULL)
  {
    printf("No X-data availablel!\n");
    exit(1);
  }

  param_validate(p);

  if(p->verbose>0)
  {
    param_show(p, stdout);
  }
  param_show(p, p->logf);

  struct sigaction act;
  memset (&act, '\0', sizeof(act));
  act.sa_handler = stoprun;
  sigaction(SIGINT, &act, NULL);

  loginfo(p, " --> Solving ... \n");

#ifdef SDL
  pthread_t th;
  liveXLview xlview;

  if(p->liveView == 1)
  {
    xlview.X = X;
    xlview.L = p->L;
    xlview.N = p->N;
    xlview.quit = 0;

    pthread_create(&th, // thread
        NULL, // pthread_attrib_t
        liveview_t, // function
        &xlview); // arg
  }
#endif

  if(p->dynamic == 1)
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
  } 
  else { 
    printf("Please use -D\n");
    exit(-1);
  }


  //radpos(X, p, stdout);

  time(&nowtime);
  p->time_final = difftime(nowtime, starttime);

  param_summary(p, X, stdout);
  param_summary(p, X, p->logf);

#ifdef SDL
  if(p->liveView == 1)
  {
    xlview.quit = 1;
    pthread_join(th, NULL);
  }
#endif

  if(p->cmmz == 1)
  {
    char * cmmfile = malloc(1024*sizeof(char));
    sprintf(cmmfile, "%s/cmmdump.cmm.gz", p->ofoldername);

    cmmwritez(cmmfile, X, p->N, p->r0, p->I, p->NI, p->L);
    free(cmmfile);
  } else {
    char * cmmfile = malloc(1024*sizeof(char));
    sprintf(cmmfile, "%s/cmmdump.cmm", p->ofoldername);

    cmmwrite(cmmfile, X, p->N, p->r0, p->I, p->NI, p->L);
    free(cmmfile);
  }
#ifndef NDEBUG
  errline(__LINE__);
#endif

  param_dumpX(p, X);

  if(1) {
    char * rpfname = malloc(1024*sizeof(char));
    sprintf(rpfname, "%s/rad.csv", p->ofoldername);
    FILE * rpfile = fopen(rpfname, "w");
    assert(rpfile != NULL);
    radpos(X, p, rpfile);
    fclose(rpfile);
    free(rpfname);
  }

  fclose(p->logf);
  param_free(p);
  free(X);

  if(0)
  {
    // Get max memory used:
    printf("sudo cat /proc/%u/status | grep VmPeak\n", getpid());
    // For 2300 points:
    // VmPeak:	   20964 kB
    sleep(30);
  }


  return 0;
}
