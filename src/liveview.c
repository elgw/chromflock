#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <SDL.h>
#include "hsvrgb.h"
#include "ellipsoid.h"

#include "liveview.h"

/* Idea:
 * Create a light weight renderer that can be plugged in anywhere to monitor the status of some
 * matrix 
 *
 *  - Use pthreads to run separately
 *  - Update window whenever the appointed data is changed
 *    or use some kind of signalling.
 *    Also add a blocking command to destroy the window (before the data is removed).
 *
 * */

float mousex = 0;
float mousey = 0;
int mousedown = 0;

static uint8_t cmap[] = {255,255,255, 	
  240,163,255,
  0,117,220, 	
  153,63,0, 	
  76,0,92, 	
  25,25,25, 	
  0,92,49, 	
  43,206,72, 	
  255,204,153, 	
  128,128,128, 	
  148,255,181, 	
  143,124,0, 	
  157,204,0, 	
  194,0,136, 	
  0,51,128, 	
  255,164,5, 	
  255,168,187, 	
  66,102,0, 	
  255,0,16, 	
  94,241,242, 	
  0,153,143, 	
  224,255,102, 	
  116,10,255, 	
  153,0,0, 	
  255,255,128, 	
  255,255,0, 	
  255,80,5};

typedef struct {
  size_t N;

  uint8_t * L; // labels
  double * X; // pointer to "live" data
  double * XZ; // [x, y, z, l]

  int window_w;
  int window_h;

  SDL_Window* window;
  SDL_Renderer* renderer;
  SDL_Surface * surface;
  SDL_Texture * texture;

  SDL_Rect * SrcR;
  SDL_Rect * DestR;

  char * title;
  int pause;

  size_t nFrames; // number of rendered frames

  int done;
  double r0;
  double * Rot; // rotation matrix
  elli * E;
  double * EA;
} scene;


typedef struct {
  SDL_Texture * texture;
} bead;

static double min_double(double a, double b)
{
  if (a < b){
    return a;
  }
  return b;
}

static int imin(int a, int b)
{
  if (a < b){
    return a;
  }
  return b;
}


static void drawBead(uint32_t * pixels, int width, int height, int label)
{

  if(width != height)
  {
    printf("Width not equal to height\n");
    exit(0);
  }

  double * RGB = malloc(3*sizeof(double));
  double * HSV = malloc(3*sizeof(double));
  double * HSV2 = malloc(3*sizeof(double));
  double * RGB2 = malloc(3*sizeof(double));

  if(label>24)
  {
    label = 0;
  }

  int64_t r = cmap[3*label];
  int64_t g = cmap[3*label+1];
  int64_t b = cmap[3*label+2];

  RGB[0] = (double) r/255.0;
  RGB[1] = (double) g/255.0;
  RGB[2] = (double) b/255.0;

  rgb2hsv(RGB, HSV);

  r = 255;
  g = 255;
  b = 0;
  int a = 0;

  int br = round((width-1)/2);

  for(int xx = -br; xx<=br; xx++)
  {
    for(int yy = -br; yy<=br; yy++)
    {
      double r = sqrt(pow(xx,2) + pow(yy, 2));

      a = 0;
      if(r <= br)
      {
        memcpy(HSV2, HSV, 3*sizeof(double));
        HSV2[2] *= sqrt((br-r)/br);
        hsv2rgb(HSV2, RGB2);

        r= round(RGB2[0]*255.0);
        g= round(RGB2[1]*255.0);
        b= round(RGB2[2]*255.0);
        a=255;
      }

      pixels[(br+xx)+width*(br+yy)] = r + 256*g + 256*256*b + 256*256*256*a;
    }
  }

  free(RGB);
  free(HSV);
  free(RGB2);
  free(HSV2);
}

void bead_init(scene * s, bead * b, int label)
{
  int width = 51;
  int height = 51;
  int depth = 32;
  int pitch = 4*width;

  uint32_t * pixels = malloc(width*height*sizeof(uint32_t));

  memset(pixels, 128, height*width*sizeof(uint32_t));
  drawBead(pixels, width, height, label);

  uint32_t rmask = 0x000000ff;
  uint32_t gmask = 0x0000ff00;
  uint32_t bmask = 0x00ff0000;
  uint32_t amask = 0xff000000;

  SDL_Surface* surf = SDL_CreateRGBSurfaceFrom(
      (void*) pixels, 
      width, 
      height, 
      depth, 
      pitch, 
      rmask, 
      gmask, 
      bmask, 
      amask);

  b->texture = SDL_CreateTextureFromSurface(
      s->renderer, 
      surf);

  SDL_FreeSurface(surf);
  free(pixels);
}

void bead_free(bead *b)
{
  SDL_DestroyTexture(b->texture);
}

static void getError()
{
  // If there is an SDL error, print it. If quit==1, set flag for
  // program termination

  if(strlen(SDL_GetError())>0)
  {
    printf("SDL ERROR\n");
    printf("%lu\n", strlen(SDL_GetError()));
    printf("%s\n", SDL_GetError());
    SDL_Delay(100);
    SDL_ClearError();
  }
}

static void gInit(scene * s)
{
  // Initialize graphics

  if (SDL_Init(SDL_INIT_VIDEO) == 0) {

    SDL_SetHint(SDL_HINT_RENDER_VSYNC, "1"); 

    // s->surface??

    s->window = SDL_CreateWindow(s->title,
        SDL_WINDOWPOS_UNDEFINED,
        SDL_WINDOWPOS_UNDEFINED,
        s->window_w, s->window_h,
        SDL_WINDOW_OPENGL + SDL_WINDOW_RESIZABLE);

    s->renderer = SDL_CreateRenderer(s->window, -1, 0);

  }
  else
  {
    getError();
    printf("Failed SDL_Init\n");
    printf("%s\n", SDL_GetError());
    SDL_Delay(10000);
  }

}


static void getEvents(scene * s)
{
  // Handle SDL events, i.e., keyboard, mouse, resize, ...


  SDL_GetWindowSize(s->window, &s->window_w, &s->window_h); // get window size

  SDL_Event evt;

  while(SDL_PollEvent(&evt)) { 

    if(evt.type == SDL_QUIT) {
      s->done = 1;
    }
    if(evt.type == SDL_MOUSEBUTTONDOWN)
    {
      mousex = (float) evt.motion.x;
      mousey = (float) evt.motion.y;
      mousedown = 1;
    }
    if(evt.type == SDL_MOUSEMOTION)
    {
      if(mousedown == 1)
      {
      float deltaX =  (float) evt.motion.x - mousex;
      float deltaY =  (float) evt.motion.y - mousey;
      mousex = (float) evt.motion.x;
      mousey = (float) evt.motion.y;
      double scale = min_double(s->window_w, s->window_h);
      rot_x(s->Rot, 4*deltaY/scale);
      rot_y(s->Rot, 4*deltaX/scale);
      }
    }

    if(evt.type == SDL_MOUSEBUTTONUP)
    {
      mousedown = 0;

      // SDL_WarpMouseInWindow(s->window, 
      //     s->window_w/2.0, s->window_h/2.0);
    }   

    if(evt.type == SDL_KEYDOWN)
    {
      if (evt.key.keysym.sym == SDLK_RIGHT || evt.key.keysym.sym == SDLK_w)
      {
        rot_y(s->Rot, 0.1);
      }
      if (evt.key.keysym.sym == SDLK_LEFT || evt.key.keysym.sym == SDLK_s)
      {
        rot_y(s->Rot, -0.1);
      }

      if (evt.key.keysym.sym == SDLK_UP || evt.key.keysym.sym == SDLK_q)
      {
        rot_x(s->Rot, 0.1);
      }
      if (evt.key.keysym.sym == SDLK_DOWN || evt.key.keysym.sym == SDLK_a)
      {
        rot_x(s->Rot, -0.1);
      }
      if (evt.key.keysym.sym == SDLK_e)
      {
        rot_z(s->Rot, 0.1);
      }
      if (evt.key.keysym.sym == SDLK_d)
      {
        rot_z(s->Rot, -0.1);
      }



      if (evt.key.keysym.sym == SDLK_ESCAPE)
      {
        s->done = 1;
      }
      if (evt.key.keysym.sym == SDLK_SPACE) {
        s->pause++;
        if(s->pause == 2)
        {
          s->pause = 0;
        }
        printf("s->space=%d\n", s->pause);
      }
    }
  }
  }

  static int zcmp(const void * A, const void * B)
  {
    double * P = (double * ) A;
    double * Q = (double * ) B;

    if(P[2] > Q[2])
      return 1;
    if(P[2] < Q[2])
      return -1;

    return 0;
  }

  static void mattrans(double * T, double * A)
  {
    T[0] = A[0]; T[3] = A[1]; T[6] = A[2];
    T[1] = A[3]; T[4] = A[4]; T[7] = A[5];
    T[2] = A[6]; T[5] = A[7]; T[8] = A[8];
  }

#if FALSE
static void matshow(double * X)
  {
    printf("[");
    for(int mm = 0; mm<3; mm++)
    {
      for(int nn = 0; nn<3; nn++)
      {
        printf("%f ", X[mm+3*nn]);
        if(nn+1 < 3)
        {
          printf(", ");
        }

      }
      printf("\n");
    }
    printf("]\n");
  }
#endif

  static void matmul(double * Z, double * X, double * Y)
  {
    // Z = X*Y, all are 3x3 matrices
    // X,Y and Z are allowed to point to the same memory
    double T[9] = {0,0,0,0,0,0,0,0,0};
    for(int mm = 0; mm<3; mm++)
    {
      for(int nn = 0; nn<3; nn++)
      {
        for(int ii = 0; ii<3; ii++)
        {
          T[mm+3*nn] += X[ii*3 + mm]*Y[nn*3 + ii];
        }
      }
    }
    memcpy(Z, T, 9*sizeof(double));
    return;
  }

  static void rot_x(double * SR, double theta)
  {
    // Rotate around x - axis
    double R[9];
    R[0] = 1; R[3] = 0;           R[6] = 0;
    R[1] = 0; R[4] = cos(theta);  R[7] = sin(theta);
    R[2] = 0; R[5] = -sin(theta); R[8] = cos(theta);
    matmul(SR, R, SR);
  }

  static void rot_y(double * SR, double theta)
  {
    // Rotate around y - axis
    double R[9];
    R[0] = cos(theta);  R[3] = 0; R[6] = sin(theta);
    R[1] = 0;           R[4] = 1; R[7] = 0;
    R[2] = -sin(theta); R[5] = 0; R[8] = cos(theta);
    matmul(SR, R, SR);
  }


  static void rot_z(double * SR, double theta)
  {
    // Rotate around z - axis
    double R[9];
    R[0] =  cos(theta); R[3] = sin(theta); R[6] = 0;
    R[1] = -sin(theta); R[4] = cos(theta); R[7] = 0;
    R[2] = 0;           R[5] = 0;          R[8] = 1;
    matmul(SR, R, SR);
  }

  static void rot_point(double * R, double * X)
  {
    // Rotate around z - axis
    double x = R[0] * X[0] + R[3]*X[1] + R[6]*X[2];
    double y = R[1] * X[0] + R[4]*X[1] + R[7]*X[2];
    double z = R[2] * X[0] + R[5]*X[1] + R[8]*X[2];
    X[0] = x;
    X[1] = y;
    X[2] = z;
  }

  static void copy_sort(scene * s)
  {
    for(size_t kk = 0; kk < s->N; kk++)
    {
      for(size_t idx = 0; idx<3; idx++)
      {
        s->XZ[4*kk + idx] = s->X[3*kk + idx];
      }
      s->XZ[4*kk + 3] = s->L[kk];
      rot_point(s->Rot, s->XZ+4*kk);
    }

    qsort(s->XZ, s->N, 4*sizeof(double), zcmp);
    return;
  }


  static void render(scene * s, bead * beads)
  {

    if(s->pause == 0)
    {
      copy_sort(s);
    }

    SDL_SetRenderDrawColor(s->renderer, 255.0, 255.0, 255.0, SDL_ALPHA_OPAQUE);
    //SDL_SetRenderDrawColor(s->renderer, 0.0, 0.0, 0.0, SDL_ALPHA_OPAQUE);
    SDL_RenderClear(s->renderer);


    // Bead radius
    double bead_radius = s->r0;
    // In terms of screen pixels
    int br = round((double) imin(s->window_w, s->window_h)*bead_radius);

    int mid = imin(s->window_w / 2, s->window_h /2);
    // Figure out offset
    int woff = 0;
    int hoff = 0l;
    int d1 = s->window_w - s->window_h;
    if(d1 > 0)      
      woff = d1/2;
    if(d1 < 0)
      hoff = -d1/2;

    /* Draw domain */
    SDL_SetRenderDrawColor(s->renderer, 0, 0, 0, 255);

    double ER[9];
    memcpy(ER, s->EA, 9*sizeof(double));
    double RT[9];
    mattrans(RT, s->Rot);
    matmul(ER, ER, RT);
    matmul(ER, s->Rot, ER);

    for(int delta = 0; delta<3; delta++)
    {
      double x0 = -10;
      double y0 = -10;

      for(double theta = 0.01; theta<=6*M_PI; theta = theta+0.02)
      {
        //      double x = mid+woff+(mid+delta)*sin(theta);
        //      double y = hoff+mid+(mid+delta)*s->E->b*cos(theta);

        double x = cos(theta); double y = sin(theta);
        double scale = // ||x^TAx||
          x*(ER[0]*x + ER[3]*y) + y*(ER[1]*x + ER[4]*y);

        x = x*(mid+delta)/sqrt(scale) + mid + woff;
        y = y*(mid+delta)/sqrt(scale) + mid + hoff;
        //      printf("Scale: %f (%f, %f)\n", scale, x, y);

        if(x0>-10)
          SDL_RenderDrawLine(s->renderer, round(x0), round(y0), round(x), round(y));
        y0 = y;
        x0 = x;
      }

      /* Draw beads */
      for(size_t kk = 0; kk<s->N; kk++)
      {
        int label = ((int) s->XZ[4*kk+3] ) % 32;
        if(label < 26)
        {
          SDL_Rect SrcR;
          SDL_Rect DestR;

          SrcR.x = 0;
          SrcR.y = 0;
          SrcR.w = 200;
          SrcR.h = 200;

          SDL_QueryTexture(beads[label].texture, NULL, NULL, &SrcR.w, &SrcR.h);


          DestR.w = br;
          DestR.h = br;
          DestR.x = mid + mid*s->XZ[4*kk] - DestR.w/2 + woff;
          DestR.y = mid + mid*s->XZ[4*kk+1] - DestR.w/2 + hoff;

          // printf("%d %d %d %d\n", DestR.x, DestR.y, DestR.w, DestR.h);

          SDL_RenderCopy(s->renderer, beads[label].texture, &SrcR, &DestR);
        }
      }

    }

    SDL_RenderPresent(s->renderer);
    s->nFrames++;
  }

  void * liveview_t(void * conf)
  {
    liveXLview * xlconf = (liveXLview * ) conf;
    liveview(xlconf->X, xlconf->L, xlconf->N, &xlconf->quit, xlconf->r0, xlconf->E);
    return NULL;
  }

  int liveview(double * X, uint8_t * L, size_t N, int * quit, double r0, elli * E)
  {

    // printf("%d\n", __LINE__); fflush(stdout);

    scene * s = malloc(sizeof(scene));

    s->pause = 0;
    s->X = X;
    s->N = N;
    s->XZ = malloc(s->N*4*sizeof(double));
    s->X = X;
    s->L = L;
    s->done = 0;
    s->window_h = 512;
    s->window_w = 512;
    s->title = malloc(512*sizeof(char));
    s->nFrames = 0;
    s->r0 = r0;
    s->E = E;
    s->Rot = calloc(9, sizeof(double));
    s->Rot[0] = 1; s->Rot[4] = 1; s->Rot[8] = 1;
    s->EA = calloc(9, sizeof(double));
    s->EA[0] = 1.0/pow(s->E->a, 2);
    s->EA[4] = 1.0/pow(s->E->b, 2);
    s->EA[8] = 1.0/pow(s->E->c, 2);
    sprintf(s->title, "live X view");

    // Initialize window
    gInit(s);

    // Initialize beads

    bead * beads = malloc(26*sizeof(bead));
    for(int bb = 0; bb<26; bb++)
    {
      bead_init(s, beads+bb, bb);
    }


    while(s->done == 0 && quit[0] == 0)
    {

      render(s, beads);
      getEvents(s);

      usleep(1000000.0/24.0);
    }

    fprintf(stdout, "Rendered %zu times\n", s->nFrames);

    SDL_DestroyRenderer(s->renderer);
    SDL_DestroyWindow(s->window);

    for(int bb = 0; bb<26; bb++)
    {
      bead_free(&beads[bb]);
    }



    free(s->XZ);
    free(s->title);
    free(s);


    return 0;
  }
