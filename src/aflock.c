/* 
 * Implementation notes:
 *
 * - Is is assumed that all coordinates will fit into memory.
 *
 * - Pairwise distances are computed when needed and note stored.
 *
 * - The contact maps are only read and saved one by one, changes to 
 *   be performed are queued and written to disk either when the queue is
 *   full or at the end.
 *
 * General:
 * - Uses infinite math so DON'T compile with -ffitnite-math-only
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include <assert.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <getopt.h>
#include <pthread.h>

#include "aflock.h"
#include "wio.h"

size_t allocated = 0;


// Activation Distance Struct
// To be passed to threads
typedef struct{
  double * A; // For storing activation distances
  double * AD; 
  fconf * fc;
  size_t thread;
  size_t nThreads;
  double th_low;
  double th_high;
  chrom * flock;
} adstruct; 



fconf * fconf_init()
{
  fconf * c = malloc(sizeof(fconf));

  c->afname = NULL;
  c->nBeads = 0;
  c->nStruct = 0;

  c->rfname = NULL;
  c->prfname = NULL;

  c->th_low = 0;
  c->th_high = 0;  
  c->mode = MODE_UNKNOWN;
  c->QS = 2500;

  c->r0 = -1;
  c->vq = -1;

  c->mflock_arguments = NULL;

  c->experimental = 0;
  c->get_quality = 0;
  int nCpus = sysconf(_SC_NPROCESSORS_ONLN);
  if(nCpus < 0)
  {
    nCpus = 4;
  }
  c->nThreads = nCpus;

  return c;
}

void fconf_free(fconf * c)
{
  if(c->A != NULL)
  {
    free(c->A);
  }
  if(c->afname != NULL)
  {
    free(c->afname);
  }
  if(c->rfname != NULL)
  {
    free(c->rfname);
  }
  if(c->prfname != NULL)
  {
    free(c->prfname);
  }
}

// 10,000 structures, 2300 beads each, 1m 32s on MBP
// 1m 35 s on erikwfractal. Requires load of ram...

static float eudist3(float * A, float * B)
{
  /* Euclidean distance between two 3D-vectors */
  return sqrt( pow(A[0]-B[0], 2) + pow(A[1]-B[1], 2) + pow(A[2]-B[2], 2));
}

static float norm3(float * X)
{
  float n = 0;
  for(size_t kk = 0; kk<3; kk++)
    n+=pow(X[kk], 2);
  return sqrt(n);
}

void fconf_load_A(fconf * fc)
{

  fprintf(stdout, "Loading A: %s ... \n", fc->afname);

  FILE * afile = fopen(fc->afname, "r");

  if(afile == NULL)
  {
    fprintf(stderr, "Failed to open file\n");
    exit(-1);
  }

  struct stat st;
  int status = stat(fc->afname, &st);
  if(status != 0)
  {
    fprintf(stderr, "Failed to get file size of %s\n", fc->afname);
    exit(-1);
  }

  size_t fsize = st.st_size;

  fc->A = malloc(fsize);
  if(fc->A == NULL)
  {
    fprintf(stderr, "Failed to allocate memory for A\n");
    exit(-1);
  }  

  size_t nRead = fread(fc->A, sizeof(double), fsize/sizeof(double), afile);

  if(nRead*sizeof(double) != fsize)
  {
    fprintf(stderr, "nRead = %zu, fsize = %zu\n", nRead, fsize);
    exit(-1);
  }

  fclose(afile);

  for(size_t kk = 0 ; kk<fsize/sizeof(double) ; kk++)
  {
    if(isfinite(fc->A[kk]))
    {
      if(fc->A[kk] < 0)
      {
        fprintf(stderr, "ERROR: A[%zu] = %f < 0\n", kk, fc->A[kk]);
        exit(-1);
      }
      if(fc->A[kk] > 1)
      {
        fprintf(stderr, "ERROR: A[%zu] = %f > 1\n", kk, fc->A[kk]);
        exit(-1);
      }
    }
  }


  fc->nBeads = sqrt(fsize/sizeof(double));
  assert(pow(fc->nBeads, 2) == fsize/sizeof(double));
  fprintf(stdout, "   .. contains [%zu x %zu] elements.\n", fc->nBeads, fc->nBeads);
  return;
}

int chrom_init(chrom * c, size_t nQ, size_t n)
{
  c->X = NULL;
  c->nQ = 0; // Number of elements IN the queue
  c->Q = malloc(2*nQ*sizeof(uint32_t));

  c->wfName = malloc(32*sizeof(char));
  sprintf(c->wfName, "cf_%06zu/W.uint8.gz", n+1);
  c->W = 0;

  if(c->wfName == NULL)
    return 1;

  c->xfName = malloc(32*sizeof(char));
  sprintf(c->xfName, "cf_%06zu/coords.csv", n+1);

  if(c->xfName == NULL)
    return 1;

  return 0;
}

void ch_free(chrom * c)
{
  if(c->W != NULL)
    free(c->W);
  if(c->X != NULL)
    free(c->X);
  if(c->Q != NULL)
    free(c->Q);
  free(c->wfName);
  free(c->xfName);
  c->nQ = 0;
  return;
}


int ch_load_X(fconf * fc, chrom * c)
{

  c->X = malloc(fc->nBeads*3*sizeof(float));
  // Read from c->xfName ...

  //  fprintf(stdout, "Reading X-data from %s\n", c->xfName);

  FILE * f = fopen(c->xfName, "r");
  if(f == NULL)
  {
    fprintf(stderr, "\rCan't open %s\n", c->xfName);
    fprintf(stderr, "Did you forget to run mflock after aflock -I?\n");
    exit(-1);
  }

  char * line = malloc(1024*sizeof(char));
  size_t len = 1024*sizeof(char);

  char delim[] = ",";
  for(size_t ll = 0; ll<fc->nBeads; ll++)
  {
    int read = getline(&line, &len, f);
    if(read == -1) 
    {
      printf("Failed to read line %zu\n", ll+1);
      exit(-1);
    }
    char *ptr = strtok(line, delim);
    c->X[3*ll] = atof(ptr);
    ptr = strtok(NULL, delim);
    c->X[3*ll+1] = atof(ptr);
    ptr = strtok(NULL, delim);
    c->X[3*ll+2] = atof(ptr);
  }

  free(line);
  fclose(f);

  return 0;
}

int cmp_float_reverse(const void * A, const void * B)
{
  float * dA = (float *) A;
  float * dB = (float *) B;

  if(dA[0] > dB[0])
    return 1;
  if(dA[0] < dB[0])
    return -1;
  return 0;
}

int cmp_float(const void * A, const void * B)
{
  float * dA = (float *) A;
  float * dB = (float *) B;

  if(dA[0] > dB[0])
    return 1;
  if(dA[0] < dB[0])
    return -1;
  return 0;
}

void chrom_assign(fconf * fc, chrom * ch, size_t kk, size_t ll)
{
  if(ch->nQ == fc->QS)
  {
    // Empty the queue first by writing to disk.
    struct_write_W(fc, ch);
  } 

  // Insert in Queue
  ch->Q[2*ch->nQ] = kk;
  ch->Q[2*ch->nQ+1] = ll;

  ch->nQ++;
}


chrom * load_structures(fconf * fc)
{
  fprintf(stdout, "Initializing %zu structures ...\n", fc->nStruct);
  // Allocate memory for their data
  chrom * flock = malloc(fc->nStruct*sizeof(chrom));

  // Load them
  for(size_t kk = 0; kk<fc->nStruct; kk++)
  {
    printf("\r%zu", kk+1); fflush(stdout);
    chrom_init(&flock[kk], fc->QS, kk);
  }
  printf("\r");

  if(fc->mode == MODE_INIT)
  {
    printf("No coordinates to be loaded\n");
  } 
  else 
  {
    printf("Loading coordinates ...\n");

    for(size_t kk = 0; kk<fc->nStruct; kk++)
    {
      printf("\r%zu", kk); fflush(stdout);
      ch_load_X(fc, &flock[kk]);
    }
    printf("\r");
  }


  printf("\r");

  fprintf(stdout, "   ... done.\n");
  return flock;
}

void writew(char * wfname, uint8_t * W, size_t nElements)
{
  wio_write(wfname, nElements, W);
}

void struct_write_W0(fconf * fc, chrom * cc, uint8_t * W0)
{
  wio_write(cc->wfName, pow(fc->nBeads,2), W0);
}

void struct_write_R0(fconf * fc, chrom * cc, double * R0)
{
  wio_write(cc->rfName, fc->nBeads*sizeof(double), (void *) R0);
}

uint8_t * readw(char * wfName, size_t nElements)
{
  return wio_read(wfName, &nElements);
}

void struct_write_W_AD(fconf * fc, 
    chrom * c, 
    double * AD, 
    double th_high, 
    double th_low)
{

  size_t nRead = pow(fc->nBeads, 2);
  uint8_t * W = readw(c->wfName, nRead);
  double * A = fc->A;

  for(size_t aa = 0; aa < fc->nBeads; aa++)
  {
    for(size_t bb = aa+1; bb < fc->nBeads; bb++)
    {
      size_t idx1 = aa+bb*fc->nBeads;
      size_t idx2 = bb+aa*fc->nBeads;
      // If within the specified theta range
      if( (A[idx1] < th_high) & (A[idx1] >= th_low) )
      {
        // And a activation distance is specified
        if(isfinite(AD[idx1]))
        {
          // Write 1 to W if the bead pair distance
          // for this structure is within the activation distance
          double d = eudist3(c->X+aa*3, c->X+bb*3);

          if(d<AD[idx1])
          {
            W[idx1] = 1;
            W[idx2] = 1;
          } else {
            W[idx1] = 0;
            W[idx2] = 0;
          }
        } else {
          W[idx1] = 0;
          W[idx2] = 0;
        }
      }
    }
  }
  // Write back to disk
  wio_write(c->wfName, nRead, W);
  free(W);
  return;
}

void struct_write_W(fconf * fc, chrom * c)
{
  //  fprintf(stdout, "Will open %s and write %zu new constraints\n", c->wfName, c->nQ);

  size_t nRead = pow(fc->nBeads, 2);
  uint8_t * W = readw(c->wfName, nRead);


  // Update W with new contacts
  for(size_t pp = 0; pp<c->nQ; pp++)
  {
    size_t kk = c->Q[2*pp];
    size_t ll = c->Q[2*pp+1];

    W[kk+fc->nBeads*ll] = 1;
    W[ll+fc->nBeads*kk] = 1;
  }

  // Write back to disk
  wio_write(c->wfName, nRead, W);
  c->nQ = 0;

  free(W);

  return;
}

void flock_reset(fconf * fc, chrom * flock, double th_high, double th_low)
{
  /* Reset all contacts in each W where th_high < A <= th_low
  */
  fprintf(stdout, "Resetting contacts where %f <= A < %f\n", th_low, th_high);
  size_t nElements = (size_t) pow(fc->nBeads, 2);

  for(size_t pp = 0; pp<fc->nStruct; pp++)
  {
    printf("\r%zu", pp); fflush(stdout);

    uint8_t * W = readw(flock[pp].wfName, nElements);

    for(size_t kk = 0; kk<nElements; kk++)
    {
      if(fc->A[kk] >= th_low && fc->A[kk] < th_high)
      {
        W[kk] = 0;
      }
    }

    writew(flock[pp].wfName, W, nElements);

    free(W);
  }
  printf("\r");
  printf("Done!  \n");


}

void flock_assign_with_w_randblock(fconf * fc, chrom * flock, size_t kk, size_t ll, float prob)
{
  /* Assign a contacts between bead kk and ll with with probability prob
  */

  size_t nChrom = fc->nStruct;
  size_t nAssign = round(prob*(float) nChrom);

  // Array of pairwise distance between point kk and ll
  float * D = malloc(nChrom*sizeof(float));
  for(size_t pp = 0; pp<nChrom; pp++)
  {
    D[pp] = eudist3(flock[pp].X+kk*3, flock[pp].X+ll*3);
    if(flock[pp].W[kk+fc->nBeads*ll] == 1)
    {
      double drand = (double) rand()/ (double) RAND_MAX;
      if(D[pp] > 2.5*fc->r0)
      {
        D[pp] = drand + 2;
      } 
      else
      {
        D[pp] = drand;
      }

    }
  }

  float * DS = malloc(nChrom*sizeof(float));
  memcpy(DS, D, nChrom*sizeof(float));

  qsort(DS, nChrom, sizeof(float), cmp_float);

  float threshold = DS[nAssign-1];
  free(DS);

  if(0)
  {
    fprintf(stdout, "Found threshold: %f\n", threshold);
  }

  size_t nAssigned = 0;

  for(size_t pp = 0; pp<nChrom; pp++)
  {
    if((D[pp]<=threshold) && (nAssigned < nAssign)) // To avoid the case where multiple structures with the same distance
    {
      flock[pp].W[kk+fc->nBeads*ll] = 1;
      flock[pp].W[ll+fc->nBeads*kk] = 1;
      nAssigned++;
    } else {
      flock[pp].W[kk+fc->nBeads*ll] = 0;
      flock[pp].W[ll+fc->nBeads*kk] = 0;
    }      
  }

  if(nAssigned != nAssign)
  {
    fprintf(stderr, "To assign: %zu\n", nAssign);
    fprintf(stderr, "Assigned: %zu\n", nAssigned);
    fprintf(stderr, "Threshold: %f\n", threshold);
    assert(0);
  }

  return;
}

void flock_assign_with_w(fconf * fc, chrom * flock, size_t kk, size_t ll, float prob)
{
  /* Assign a contacts between bead kk and ll with with probability prob
  */

  size_t nChrom = fc->nStruct;
  size_t nAssign = round(prob*(float) nChrom);

  // Array of pairwise distance between point kk and ll
  float * D = malloc(nChrom*sizeof(float));
  for(size_t pp = 0; pp<nChrom; pp++)
  {
    D[pp] = eudist3(flock[pp].X+kk*3, flock[pp].X+ll*3);
    if(flock[pp].W[kk+fc->nBeads*ll] == 1)
    {
      if(D[pp] > 2.5*fc->r0)
      {
        D[pp] += 2;
      }
    }
  }

  float * DS = malloc(nChrom*sizeof(float));
  memcpy(DS, D, nChrom*sizeof(float));

  qsort(DS, nChrom, sizeof(float), cmp_float);

  float threshold = DS[nAssign-1];
  free(DS);

  if(0)
  {
    fprintf(stdout, "Found threshold: %f\n", threshold);
  }

  size_t nAssigned = 0;

  for(size_t pp = 0; pp<nChrom; pp++)
  {
    if((D[pp]<=threshold) && (nAssigned < nAssign)) // To avoid the case where multiple structures with the same distance
    {
      flock[pp].W[kk+fc->nBeads*ll] = 1;
      flock[pp].W[ll+fc->nBeads*kk] = 1;
      nAssigned++;
    } else {
      flock[pp].W[kk+fc->nBeads*ll] = 0;
      flock[pp].W[ll+fc->nBeads*kk] = 0;
    }      
  }

  if(nAssigned != nAssign)
  {
    fprintf(stderr, "To assign: %zu\n", nAssign);
    fprintf(stderr, "Assigned: %zu\n", nAssigned);
    fprintf(stderr, "Threshold: %f\n", threshold);
    assert(0);
  }

  return;
}

void flock_assign(fconf * fc, chrom * flock, size_t kk, size_t ll, float prob)
{
  /* Assign a contacts between bead kk and ll with with probability prob
  */

  size_t nChrom = fc->nStruct;
  size_t nAssign = round(prob*(float) nChrom);
  assert(nAssign > 0);
  assert(nAssign <= nChrom);
  if(nAssign == 0)
  {
    return;
  }

  if(0)
  {
    printf("prob=%f -> %zu/%zu structures\n", prob, nAssign, nChrom);
  }
  float * D = malloc(nChrom*sizeof(float));
  for(size_t pp = 0; pp<nChrom; pp++)
  {
    D[pp] = eudist3(flock[pp].X+kk*3, flock[pp].X+ll*3);
    // If we knew which contacts that failed, we could increase the distance of them here in order not to 
    // try them again... (at least not in the next round).
  }

  float * DS = malloc(nChrom*sizeof(float));
  memcpy(DS, D, nChrom*sizeof(float));

  qsort(DS, nChrom, sizeof(float), cmp_float);


  float threshold = DS[nAssign-1];

  free(DS);

  if(0)
  {
    fprintf(stdout, "Found threshold: %f\n", threshold);
  }

  size_t nAssigned = 0;

  /* (1) In case that the distances are equal for several structures
   * the contact is assigned to the first of the structures
   */

  if(fc->experimental == 0)
  {
    for(size_t pp = 0; pp<nChrom; pp++)
    {
      if((D[pp] <= threshold) &&
          (nAssigned < nAssign)) // (1)
      {
        chrom_assign(fc, flock+pp, kk, ll);
        nAssigned++;
      }
    }
  } else {
    // Experimental random assignment!

    uint8_t * set = malloc(nChrom*sizeof(uint8_t));
    memset(set, 0, nChrom*sizeof(uint8_t));

    for(size_t pp =0; pp<nAssign; pp++)
    {
      int isSet = 0;
      while(isSet == 0)
      {
        size_t pos = rand() % nChrom;
        if(set[pos] == 0)
        {
          isSet = 1;
          set[pos] = 1;
          chrom_assign(fc, flock+pos, kk, ll);
          nAssigned++;
        }
      }
    }
  }
  assert(nAssigned == nAssign);
  free(D);
  return;
}


static void usage()
{
  printf("Accepted arguments:\n");
  printf(" -A file.dat, --Afile file.dat\n\tcontact probability matrix, square double.\n");
  printf(" -n nStruct, --nStruct nStruct\n\tspecify the number of structures\n");
  printf(" -h th, --high th\n\thigher threshold\n");
  printf(" -l th, --low th\n\tlower threshold\n");
  printf(" -R radius, --radius radius\n\tgive bead radius\n");
  printf(" -Q vq, --vq vq\n\tgive volume quotient of beads/nuclei (instead of -R)\n");
  printf(" -P args, --mArgs args\n\tcommand line arguments to pass to mflock\n");
  printf(" --rpos RPOS.double\n\tsupply radial preference for each bead\n");
  printf(" --prpos PRPOS.double\n\tsupply the probability that each bead is given the radial constraint\n");
  printf("MODES:\n");
  printf(" -I, --init\n\tinitialize new structures\n");
  printf(" -F, --final\n\tfinal stuff\n");
  printf(" -U, --update\n\tupdate\n");
  printf(" -B, --blockupdate\n\tupdate and block failed contacts\n");
  printf(" -X, --update_experimental\n\texperimental update algorithm\n");
  printf(" -S, --score\n\tscore current quality (not for -I)\n");

}

static int argparsing(fconf * p, int argc, char ** argv)
{

  struct option longopts[] = {
    // MODES
    { "help",         no_argument,       NULL,   'u' },
    { "init",         no_argument,       NULL,   'I' },
    { "update",       no_argument,       NULL,   'U' },
    { "blockupdate",       no_argument,       NULL,   'B' },
    { "final",        no_argument,       NULL,   'F' },
    { "score",        no_argument,       NULL,   'S' },
    { "update_experimental", no_argument, NULL,  'X' },
    // Data
    { "Afile",        required_argument, NULL,   'A' },
    { "Rfile",        required_argument, NULL,   'R' },
    // Settings
    { "mArgs",        required_argument, NULL,   'P' },
    { "nStruct",      required_argument, NULL,   'n' },
    { "high",         required_argument, NULL,   'h' },
    { "low",          required_argument, NULL,   'l' },
    { "radius",       required_argument, NULL,   'R' },
    { "vq",           required_argument, NULL,   'Q' },
    {"rpos",          required_argument, NULL,   'r' },
    {"prpos",          required_argument, NULL,  'p' },
    {"version",     no_argument,       NULL,   'i' }, 
    { NULL,           0,                 NULL,   0   }    
  };

  int ch;
  while((ch = getopt_long(argc, argv, "A:FIh:l:n:R:P:Q:ESUBXur:p:", longopts, NULL)) != -1)
  {
    switch(ch) {
      case 'i':
        printf("Build date: %s, %s\n", __DATE__, __TIME__);
        printf("GIT HASH: %s\n", GIT_VERSION);
        printf("Compiler: %s\n", CC_VERSION);
        exit(0);
      case 'A': // sphere contact probability matrix
        p->afname = malloc(strlen(optarg)+1);
        strcpy(p->afname, optarg);
        break;
      case 'I': // init, i.e., create W
        p->mode = MODE_INIT;
        break;
      case 'n':
        p->nStruct = atol(optarg);
        break;
      case 'u':
        usage();
        exit(0);
      case 'h': // high threshold
        p->th_high = atof(optarg);
        break;
      case 'l': // low threshold
        p->th_low = atof(optarg);
        break;
      case 'R':
        p->r0 = atof(optarg);
        break;
      case 'F':
        p->mode = MODE_FINAL;
        break;
      case 'U':
        p->mode = MODE_UPDATE;
        break;
      case 'B':
        p->mode = MODE_UPDATE_BLOCKED;
        break;
      case 'X':
        p->mode = MODE_EXPERIMENTAL;
        break;
      case 'P':
        p->mflock_arguments = malloc(strlen(optarg)+1);
        strcpy(p->mflock_arguments, optarg);
        break;
      case 'r':
        p->rfname = malloc(strlen(optarg)+1);
        strcpy(p->rfname, optarg);
        break;
      case 'p':
        p->prfname = malloc(strlen(optarg)+1);
        strcpy(p->prfname, optarg);
        break;
      case 'Q':
        p->vq = atof(optarg);
        break;
      case 'E':
        p->experimental = 1;
        break;
      case 'S':
        p->get_quality = 1;
        break;
      default:
        return(1);
    }
  }

  if(p->afname == NULL)
  {
    printf("No A file was specified\n");
    return 1;
  }

  if(p->nStruct == 0)
  {
    printf("Don't know how many structures to deal with\n");
    return 1;
  }

  if((p->rfname != NULL) & (p->prfname == NULL))
  {
    printf("--rpos given but not --prpos\n");
    return 1;
  }

  if((p->rfname == NULL) & (p->prfname != NULL))
  {
    printf("--prpos given but not --rpos\n");
    return 1;
  }

  return 0;
}

void struct_add_radius(fconf * fc, chrom * cc, double * R)
{
  for(size_t pp = 0; pp<fc->nBeads; pp++)
  {
    R[pp]+=norm3(cc->X+3*pp);
  }
}

void struct_get_interactions(fconf * fc, chrom * cc, double * M, 
    double cDistance) // Contact distance in normalized units
{
  size_t N = fc->nBeads;

  for(size_t kk = 0; kk<N; kk++)
  {
    for(size_t ll = kk+1; ll<N; ll++)
    {
      if( eudist3(cc->X + 3*kk, cc->X + 3*ll) <= cDistance )
      {
        M[kk + ll*N]++;
        M[ll + kk*N]++;
      }
    }
  }
  return;
}


void set_bead_radius(fconf * fc)
{

  if((fc->vq > 0) && (fc->r0 > 0))
  {
    fprintf(stderr, "Either -R or -Q has to be specified, not both!\n");
    exit(-1);
  }

  double captureFactor = 4;

  if(fc->vq > 0)
  {
    printf("Setting beads radius from volume quotient\n");

    if( (fc->vq <= 0) || (fc->vq > 1))
    {
      fprintf(stderr, "The volume quotient has be be in ]0, 1], use .4 if unsure\n");
      exit(-1);
    }

    double vSphere = 4.0/3.0*M_PI;
    double r0 = cbrt(fc->vq/(double) fc->nBeads);
    double vBeads = (double) fc->nBeads*pow(r0, 3)*4.0/3.0*M_PI;
    fc->r0 = r0;
    fc->dContact = captureFactor*fc->r0;
    fprintf(stdout, "   Bead radius: %f\n", fc->r0);
    fprintf(stdout, "   Total bead volume: %f\n", vBeads);
    fprintf(stdout, "   Sphere volume: %f\n", vSphere);
    fprintf(stdout, "   Volume Quotient: %f\n", vBeads/vSphere);
    fprintf(stdout, "   Contact distance (%.1f x bead radius): %f\n", captureFactor, fc->dContact);
    return;
  }

  if(fc->r0 > 0)
  {
    printf("Setting volume quotient from bead radius\n");
    if(fc->r0 > 1)
    {
      fprintf(stderr, "A radius of %f does not make sense\n", fc->r0);
      exit(-1);
    }

    double vSphere = 3.0/4.0*M_PI;
    double r0 = fc->r0;
    double vBeads = fc->nBeads*pow(r0, 3)*3.0/4.0*M_PI;
    fc->vq = vBeads/vSphere;
    fc->r0 = r0;
    fc->dContact = captureFactor*fc->r0;
    fprintf(stdout, "   Bead radius: %f\n", fc->r0);
    fprintf(stdout, "   Total bead volume: %f\n", vBeads);
    fprintf(stdout, "   Sphere volume: %f\n", vSphere);
    fprintf(stdout, "   Volume Quotient: %f\n", vBeads/vSphere);
    fprintf(stdout, "   Contact distance (%.1f x bead radius): %f\n", captureFactor, fc->dContact);
    return;
  }

  fprintf(stderr, "Either -R or -Q has to be specified\n");
  exit(-1);
}

double structQuality(fconf * fc, chrom * cc)
{

  size_t W0S0 = 0;
  size_t W0S1 = 0;
  size_t W1S0 = 0;
  size_t W1S1 = 0;

  // load W
  uint8_t * W = readw(cc->wfName, pow(fc->nBeads, 2));

  size_t N = fc->nBeads;

  for(size_t kk = 0; kk<N; kk++)
  {
    for(size_t ll = kk+1; ll<N; ll++)
    {
      int contact = 0; // if the beads are in contact
      int set = 0; // if there should be a contact

      double d = eudist3(cc->X + 3*kk, cc->X + 3*ll);

      if(d < fc->dContact)
      {
        contact = 1;
      }

      if(W[kk*N+ll] == 1)
      {
        set = 1;
      }

      if(set == 1)
      { 
        if(contact == 1)
        {
          W1S1++;
        } else {
          W1S0++;
        }
      } else {
        if(contact == 1)
        {
          W0S1++;
        } else 
        {
          W0S0++;
        }
      }
    }

  }

  free(W);

  printf("W==0,S==0, %zu, W==0,S==1, %zu, W==1, S==0, %zu, W==1, S==1 %zu\n", W0S0, W0S1, W1S0, W1S1);

  //return ((double) W0S0 + (double) W1S1)/((double) W0S0 + (double) W0S1 + (double) W1S0 + (double) W1S1);
  return (double) W1S1 / ( (double) W1S0 + (double) W1S1 ); // got / wanted
}

double getQuality(fconf * fc, chrom * flock)
{

  double * quality = malloc(fc->nStruct*sizeof(double));
  memset(quality, 0, fc->nStruct*sizeof(double));

  double q = 0;

  for(size_t cc = 0; cc < fc->nStruct; cc++)
  {
    quality[cc] = structQuality(fc, flock+cc);
    q+=quality[cc];
  }

  free(quality);
  return q / fc->nStruct;
}

void echo_args(int argc, char ** argv)
{
  for(int kk = 0; kk<argc; kk++)
  {
    printf("%s ", argv[kk]);
  }
  printf("\n");
  return;
}

void calc_activation_distance_th( adstruct * s)
{
  fconf * fc = s->fc;
  double * A = s->A;
  double * AD = s->AD;
  size_t thread = s->thread;
  size_t nThreads = s->nThreads;
  double th_low = s->th_low;
  double th_high = s->th_high;
  chrom * flock = s->flock;

  float * DS = malloc(fc->nStruct*sizeof(float)); // temporary array for each bead pair
  for(size_t aa = thread; aa < fc->nBeads; aa = aa+nThreads)
  {

    for(size_t bb = aa+1; bb < fc->nBeads; bb++)
    {
      size_t idx1 = aa+bb*fc->nBeads;
      size_t idx2 = bb+aa*fc->nBeads;

      if((A[idx1] < th_high) & (A[idx1] >= th_low))
      {

        size_t nAssign = round(A[idx1]*(float) fc->nStruct);

        if(nAssign>0)
        {
          // Get distance between bead aa and bb in all structures
          for(size_t ss = 0; ss< fc->nStruct; ss++)
          {
            DS[ss] = eudist3(flock[ss].X+aa*3, flock[ss].X+bb*3);
          }
          qsort(DS, fc->nStruct, sizeof(float), cmp_float);
          double act_dist = 10;
          if(nAssign < fc->nStruct)
          {
            act_dist = (DS[nAssign-1] + DS[nAssign])/2;
          }
          AD[idx1] = act_dist;
          AD[idx2] = act_dist;
        } 
      }
    }
  }

  free(DS);
}

void * run_calc_activation_distance_th(void * P)
{
  fflush(stdout);
  adstruct * s = (adstruct *) P;
  calc_activation_distance_th(s);
  return NULL;
}


void calc_activation_distance(fconf * fc, chrom * flock, double th_high, double th_low, double * AD)
{
  /* Consider all positions in PR where p<th_high and p>= th_low
   * For each beads:
   * hand the constrains to the structures that has the closest radii already
   * similar to how the contact restraints are handled
   */

  printf("%s\n", __func__);
  double * A = fc->A;

  pthread_t * threads = malloc(fc->nThreads*sizeof(pthread_t));
  adstruct ** S = malloc(fc->nThreads*sizeof(adstruct * ));
  for(size_t kk = 0; kk<fc->nThreads; kk++)
  {
    S[kk] = malloc(sizeof(adstruct));
    S[kk]->thread = kk;
    S[kk]->nThreads = fc->nThreads;
    S[kk]->fc = fc;

    S[kk]->AD = AD;
    S[kk]->A = A;
    S[kk] -> th_high = th_high;
    S[kk] -> th_low = th_low;
    S[kk]->flock = flock;
    pthread_create(&threads[kk], 
        NULL, 
        run_calc_activation_distance_th, 
        (void *) S[kk]);
  }

  for(size_t kk = 0; kk<fc->nThreads; kk++)
  {
    pthread_join(threads[kk], NULL);
    free(S[kk]);
  }
  free(S);
  free(threads);

  fprintf(stdout, "\n\n");
  fprintf(stdout, "Writing activation distances to 'activationDistance.double'\n");
  wio_write("activationDistance.double", 
      fc->nBeads*fc->nBeads*sizeof(double),
      (void *) AD);
  fprintf(stdout, "Done\n"); fflush(stdout);

  return;
}

void flock_updateR(fconf * fc, chrom * flock, double th_high, double th_low)
{
  /* Consider all positions in PR where p<th_high and p>= th_low
   * For each beads:
   * hand the constrains to the structures that has the closest radii already
   * similar to how the contact restraints are handled
   */

  // To start with, all coordinates should already be available in flock.

  // Load preferred radial positions R
  // Load probability/proportion of structures having the constraint into PR

  printf("%s\n", __func__);

  // Read R: radius per bead and PR: probability of use of each bead
  size_t nbytes_r = 0;
  size_t nbytes_pr = 0;
  double * R = (double *) wio_read(fc->rfname, &nbytes_r);
  double * PR = (double *) wio_read(fc->prfname, &nbytes_pr);

  assert(nbytes_r == nbytes_pr);
  assert(nbytes_r == fc->nBeads*sizeof(double));

  // Find a threshold for each positions
  double * TH = (double * ) malloc(fc->nBeads*sizeof(double));
  float * DS = malloc(fc->nStruct*sizeof(float));
  for(size_t pp = 0; pp < fc->nBeads; pp++)
  {
    if((PR[pp] < th_high) & (PR[pp] >= th_low))
    {
      size_t nAssign = round(PR[pp]*(float) fc->nStruct);
      if(nAssign>0)
      {

        for(size_t ss = 0; ss< fc->nStruct; ss++)
        {
          DS[ss] = norm3(flock[ss].X+pp*3);
          //printf("%f ", DS[ss]);
        }
        //printf("\n");
        qsort(DS, fc->nStruct, sizeof(float), cmp_float);
        TH[pp] = DS[nAssign-1];
        if(0){
          printf("PR[%zu]=%f, R[%zu]=%f. nassign: %zu. Th[%zu]=%f\n", 
              pp, PR[pp], pp, R[pp], nAssign, pp, TH[pp]);
        }
      }
    }
  }
  // Read/Write all individual R-files.
  for(size_t ss = 0; ss < fc->nStruct ; ss++)
  {
    flock[ss].rfName = malloc(1024*sizeof(char));
    sprintf(flock[ss].rfName, "cf_%06zu/radius.double.gz", ss+1);

    size_t bytesRead = 0;
    double * SR = (double *) wio_read( flock[ss].rfName, &bytesRead );
    if(bytesRead != sizeof(double)*fc->nBeads)
    {
      fprintf(stderr, "%d ERROR: Read %zu bytes from %s, expected %zu.\n", __LINE__, 
          bytesRead, flock[ss].rfName, fc->nBeads*sizeof(double));
      exit(-1);
    }

    // Go through all pp "beads" of structure ss
    for(size_t pp = 0 ; pp < fc->nBeads ; pp++)
    {
      if((PR[pp] < th_high) & (PR[pp] >= th_low)) // Don't change the values for the other beads
      {
        double r = norm3(flock[ss].X+pp*3);
        if(r <= TH[pp])
        {
          SR[pp] = R[pp];
        } else {
          SR[pp] = NAN;
        }
      }   
    }
    wio_write(flock[ss].rfName, 
        fc->nBeads*sizeof(double), 
        (void *) SR);
    free(SR);
    free(flock[ss].rfName);
  }

  free(DS);
  free(TH);
  free(PR);
  free(R);
  return;
}

int main(int argc, char ** argv)
{

  echo_args(argc, argv);

  fconf * fc = fconf_init();

  int ok = argparsing(fc, argc, argv);

  if(fc->mode == MODE_UNKNOWN)
  {
    printf("No MODE specified!\n");
    exit(1);
  }

  if(ok != 0)
  {
    usage();
    exit(1);
  }

  /* Load the contact probability map specified with -A
  */

  fconf_load_A(fc);

  // set r0(qc) or qc(r0)
  // Can only be done when it is known how many beads
  // that are in the structures
  set_bead_radius(fc);

  chrom * flock = load_structures(fc);

  if(fc->get_quality == 1)
  {
    double quality = getQuality(fc, flock);
    printf("Quality: %f\n", quality);
  }

  if(fc->mode == MODE_INIT)
  {
    /* Creates initial W matrices for each structure
     * using only the contacts with theta == 1
     */
    fprintf(stdout, "Initialization mode.\n");
    fprintf(stdout, "Creating `mflock_jobs' to run with GNU parallel\n");
    FILE * jobFile = fopen("mflock_jobs", "w");
    if(jobFile == NULL)
    {
      fprintf(stderr, "Failed to create 'mflock_jobs'\n");
      exit(1);
    }

    fprintf(stdout, "Creating initial contact indication matrices, W\n");
    char * dir = malloc(32*sizeof(char));

    for(size_t kk = 0; kk<fc->nStruct; kk++)
    {
      sprintf(dir, "cf_%06zu/", kk+1);

      if(fc->rfname == NULL)
      {
        fprintf(jobFile, "mflock -w %sW.uint8.gz -o %s -x %scoords.csv -R %f", 
            dir, dir, dir, fc->r0);
      } else {
        fprintf(jobFile, "mflock -w %sW.uint8.gz -o %s -x %scoords.csv -r %sradius.double.gz -R %f", 
            dir, dir, dir, dir, fc->r0);
      }

      if(fc->mflock_arguments != NULL)
      {
        fprintf(jobFile, " %s", fc->mflock_arguments);
      }

      fprintf(jobFile, "\n");

      struct stat st;
      memset(&st, 0, sizeof(stat));

      if (stat(dir, &st) == -1) 
      {
        mkdir(dir, 0700);
      }
    }
    free(dir);
    fclose(jobFile);

    // Create the common start contacts
    uint8_t * W0 = malloc(pow(fc->nBeads, 2)*sizeof(uint8_t));
    for(size_t kk = 0; kk<pow(fc->nBeads,2); kk++)
    {
      if(fc->A[kk] == 1)
      {
        W0[kk] = 1;
      } else {
        W0[kk] = 0;
      }
    }

    for(size_t pp = 0; pp < fc->nStruct; pp++)
    {
      printf("\r%zu", pp); fflush(stdout);
      struct_write_W0(fc, &flock[pp], W0);
    }
    printf("\r");

    // Create initial radial constraints
    if(fc->rfname != NULL)
    {      
      size_t nel_r = 0;
      size_t nel_pr = 0;
      double * R = (double * ) wio_read(fc->rfname, &nel_r);
      double * PR = (double * ) wio_read(fc->prfname, &nel_pr);
      if(nel_r != nel_pr)
      {
        fprintf(stderr, "ERROR: Number of elements in %s and %s does not match\n", fc->rfname, fc->prfname);
        exit(-1);
      }
      if(nel_r != fc->nBeads*sizeof(double))
      {
        fprintf(stderr, "Error: Number of elements in %s does not match the number of beads\n", fc->rfname);
        exit(-1);
      }

      size_t used_R = 0; 
      double * R0 = malloc(fc->nBeads*sizeof(double));
      for(size_t kk = 0; kk< fc->nBeads; kk++)
      {
        R0[kk] = NAN;
        if(PR[kk] == 1)
        {
          R0[kk] = R[kk];
          used_R ++;
        }
      }
      printf("Using %zu elements of %s initially\n", used_R, fc->rfname);

      for(size_t pp = 0; pp < fc->nStruct; pp++)
      {
        printf("\r%zu", pp); fflush(stdout);
        sprintf(dir, "cf_%06zu/", pp+1);
        flock[pp].rfName = malloc(128*sizeof(char));
        sprintf(flock[pp].rfName, "%sradius.double.gz", dir);
        struct_write_R0(fc, &flock[pp], R0);
        free(flock[pp].rfName);
        flock[pp].rfName = NULL;
      }

      printf("\r");
      free(R0);
    }

    free(W0);
    printf("Done!\n");
  }


  if(fc->mode == MODE_EXPERIMENTAL)
  {
    fprintf(stdout, "Update/assignment mode. Experimental/Blocking/Random\n");
    /*
     * Removes failed contacts and assigns them to random structure
     */

    double th_high = fc->th_high;
    double th_low = fc->th_low;

    size_t nElements = fc->nBeads*fc->nBeads;

    for(size_t pp = 0; pp < fc->nStruct; pp++)
    {
      flock[pp].W = readw(flock[pp].wfName, nElements);
    }

    size_t nA = 0;
    for(size_t ll = 0; ll<fc->nBeads; ll++)
    {
      for(size_t kk = ll+1; kk<fc->nBeads; kk++)
      {
        double prob = fc->A[kk+fc->nBeads*ll];

        if( (prob >= th_low) && (prob < th_high) )
        {
          /* Assign contact between beads kk and ll
           * with probability prob to the structures in flock */
          flock_assign_with_w_randblock(fc, flock, kk, ll, prob);
          nA++;
        }
      }
    }

    fprintf(stdout, "Writing to disk\n");
    for(size_t pp = 0; pp < fc->nStruct; pp++)
    {
      printf("\r%zu", pp); fflush(stdout);
      struct_write_W0(fc, &flock[pp], flock[pp].W);
    }
    printf("\r");
    printf("\n");
  }

  if(fc->mode == MODE_UPDATE_BLOCKED)
  {
    fprintf(stdout, "Update/assignment mode. (Block failed contacts)\n");
    /* Adds new contact to W of each structure.
     * Note that there is no check that the contact
     * within the specified theta range isn't already in use, i.e.,
     * there will be confusion if aflock is run more than once for 
     * any interval of theta.
     */

    double th_high = fc->th_high;
    double th_low = fc->th_low;

    size_t nElements = fc->nBeads*fc->nBeads;

    for(size_t pp = 0; pp < fc->nStruct; pp++)
    {
      flock[pp].W = readw(flock[pp].wfName, nElements);
    }

    size_t nA = 0;
    for(size_t ll = 0; ll<fc->nBeads; ll++)
    {
      for(size_t kk = ll+1; kk<fc->nBeads; kk++)
      {
        double prob = fc->A[kk+fc->nBeads*ll];

        if( (prob >= th_low) && (prob < th_high) )
        {
          /* Assign contact between beads kk and ll
           * with probability prob to the structures in flock */
          flock_assign_with_w(fc, flock, kk, ll, prob);
          nA++;
        }
      }
    }

    fprintf(stdout, "Writing to disk\n");
    for(size_t pp = 0; pp < fc->nStruct; pp++)
    {
      printf("\r%zu", pp); fflush(stdout);
      struct_write_W0(fc, &flock[pp], flock[pp].W);
    }
    printf("\r");
    printf("\n");
  }

  if(fc->mode == MODE_UPDATE)
  {
    fprintf(stdout, "Update/assignment mode.\n");
    /* Adds new contact to W of each structure.
    */

    double th_high = fc->th_high;
    double th_low = fc->th_low;

    if(fc->rfname != NULL)
    {
      flock_updateR(fc, flock, th_high, th_low);
    }

    // Activation distance initialize as NAN
    double * AD = malloc(fc->nBeads*fc->nBeads*sizeof(double));
    for(size_t kk = 0 ; kk<pow(fc->nBeads,2) ; kk++)
    {
      AD[kk] = NAN;
    }

    fprintf(stdout, "Calculating the activation distances\n");
    calc_activation_distance(fc, flock, th_high, th_low, AD);

    fprintf(stdout, "Writing to disk\n");
    for(size_t kk = 0; kk<fc->nStruct; kk++)
    {
      if(kk%10 == 0)
      {
        printf("\r%zu", kk); fflush(stdout);
      }
      struct_write_W_AD(fc, &flock[kk], AD, th_high, th_low);
      ch_free(&flock[kk]);
    }
    printf("\r");
    printf("\n");
  }

  if(fc->mode == 999) // OLD MODE_UPDATE
  {
    fprintf(stdout, "Update/assignment mode.\n");
    /* Adds new contact to W of each structure.
     * Note that there is no check that the contact
     * within the specified theta range isn't already in use, i.e.,
     * there will be confusion if aflock is run more than once for 
     * any interval of theta.
     */

    double th_high = fc->th_high;
    double th_low = fc->th_low;

    flock_reset(fc, flock, th_high, th_low);

    if(fc->rfname != NULL)
    {
      flock_updateR(fc, flock, th_high, th_low);
    }

    size_t nA = 0;
    size_t mpos = 0;
    for(size_t ll = 0; ll<fc->nBeads; ll++)
    {
      for(size_t kk = ll+1; kk<fc->nBeads; kk++)
      {
        mpos++;
        if(mpos%1000 == 0)
        {
          printf("\r %zu/%zu", mpos, fc->nBeads*(fc->nBeads-1)/2);
        }

        double prob = fc->A[kk+fc->nBeads*ll];

        if( (prob >= th_low) && (prob < th_high) )
        {
          /* Assign contact between beads kk and ll
           * with probability prob to the structures in flock */
          flock_assign(fc, flock, kk, ll, prob);
          nA++;
        }
      }
    }
    printf("\r\n");

    printf("Found %zu contact pairs with p in ]%f, %f] \n", nA, th_high, th_low);

    fprintf(stdout, "Writing to disk\n");
    for(size_t kk = 0; kk<fc->nStruct; kk++)
    {
      printf("\r%zu", kk); fflush(stdout);
      struct_write_W(fc, &flock[kk]);
      ch_free(&flock[kk]);
    }
    printf("\r");
    printf("\n");
  }

  if(fc->mode == MODE_FINAL)
  {
    /* Loops over all the structures and 
     * creates a joint contact map.
     * Note that the capture distance, cDistance is calculated 
     * as a function of the bead radius. The bead radius
     * in this step does not have to match the one used in the simulations.
     */

    fprintf(stdout, "Finalization mode.\n");

    fprintf(stdout, "Capturing contacts, where d<%f\n", fc->dContact);
    // Generate contact probability matrix etc

    size_t msize = pow(fc->nBeads, 2)*sizeof(double);
    double * rprof = malloc(fc->nBeads*sizeof(double));
    memset(rprof, 0, fc->nBeads*sizeof(double));

    double * M = malloc(msize);
    memset(M, 0, msize);

    for(size_t pp = 0; pp < fc->nStruct; pp++)
    {
      printf("\r%zu",pp+1);
      struct_get_interactions(fc, &flock[pp], M, fc->dContact);
      struct_add_radius(fc, &flock[pp], rprof);
      fflush(stdout);
    }

    printf("\r");
    fflush(stdout);

    char * mfilename = malloc(1024*sizeof(char));
    sprintf(mfilename, "all_contacts.double");
    fprintf(stdout, "Writing contact map to: %s\n", mfilename);
    FILE * mout = fopen(mfilename, "w");

    if(mout == NULL)
    {
      printf("Failed to open output file %s\n", mfilename);
      exit(1);
    }

    size_t nwrite = fwrite(M, sizeof(double), (size_t) pow(fc->nBeads,2), mout);
    printf("Wrote %zu doubles \n", nwrite);
    assert(nwrite == (size_t) pow(fc->nBeads,2));

    char * rproffname = malloc(1024*sizeof(char));
    sprintf(rproffname, "radial_profile.csv");
    fprintf(stdout, "%s\n", rproffname);
    FILE * rproffile = fopen(rproffname, "w");
    for(size_t pp = 0; pp < fc->nBeads; pp++)
    {
      fprintf(rproffile, "%f\n", rprof[pp]/(double) fc->nStruct);
    }
    free(rprof);

    printf("Summing all W into Wtotal.double\n");
    double * Wtotal = malloc(pow(fc->nBeads,2)*sizeof(double));
    memset(Wtotal, 0, pow(fc->nBeads,2)*sizeof(double));
    char * wFileName = malloc(1024*sizeof(char));
    for(size_t pp = 0; pp < fc->nStruct ; pp++)
    {
      sprintf(wFileName, "cf_%06zu/W.uint8.gz", pp+1);
      printf("\r%s", wFileName);
      size_t nRead = 0;
      uint8_t * W = (uint8_t *) wio_read(wFileName, &nRead);

      for(size_t qq = 0; qq< pow(fc->nBeads, 2); qq++)
      {
        Wtotal[qq]+=(double) W[qq];
      }
      free(W);
    }
    free(wFileName);
    printf("\n\n");

    fprintf(stdout, "Writing\n");
    wio_write("Wtotal.double", pow(fc->nBeads,2)*sizeof(double), (uint8_t *) Wtotal);
    free(Wtotal);

    fclose(rproffile);
    free(rproffname);

    fclose(mout);
    free(mfilename);

    free(M);
    for(size_t kk = 0; kk<fc->nStruct; kk++)
    {
      ch_free(&flock[kk]);
    }

  }


  free(flock);

  fconf_free(fc);
  free(fc);

  return 0;
}
