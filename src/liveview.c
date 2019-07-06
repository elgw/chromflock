#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>
#include <SDL.h>
#include "hsvrgb.h"

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
} scene;


typedef struct {
  SDL_Texture * texture;
} bead;

static int imin(int a, int b)
{
  if (a < b)
    return a;
  return b;
}

static void drawBead(uint32_t * pixels, int width, int height, int label)
{


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

void bead_free(scene * s, bead *b)
{
  SDL_DestroyTexture(b->texture);
}

static void getError(scene * s)
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
    getError(s);
    printf("Failed SDL_Init\n");
    printf("%s\n", SDL_GetError());
    SDL_Delay(10000);
  }

  getError(s);
}


static void getEvents(scene * s)
{
  // Handle SDL events, i.e., keyboard, mouse, resize, ...


  SDL_GetWindowSize(s->window, &s->window_w, &s->window_h); // get window size

  SDL_Event evt;

  //while(SDL_PollEvent(&evt)) { // TODO: has to be in main thread
  while(0) {

    if(evt.type == SDL_QUIT) {
      s->done = 1;
    }

    if(evt.type == SDL_KEYDOWN)
    {
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
  if(P[2] > Q[2])
    return -1;

  return 0;
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


  int br = round(imin(s->window_w, s->window_h)*0.040356/2.0);
  int mid = imin(s->window_w / 2, s->window_h /2);
  
  for(size_t kk = 0; kk<s->N; kk++)
  {
    int label = s->XZ[4*kk+3];
  if(label < 26)
  {
    SDL_Rect SrcR;
    SDL_Rect DestR;

    SrcR.x = 0;
    SrcR.y = 0;
    SrcR.w = 200;
    SrcR.h = 200;

    SDL_QueryTexture(beads[label].texture, NULL, NULL, &SrcR.w, &SrcR.h);


    DestR.w = 2*br+1;
    DestR.h = 2*br+1;
    DestR.x = mid + mid*s->XZ[4*kk] - DestR.w/2;
    DestR.y = mid + mid*s->XZ[4*kk+1] - DestR.w/2;

   // printf("%d %d %d %d\n", DestR.x, DestR.y, DestR.w, DestR.h);

    SDL_RenderCopy(s->renderer, beads[label].texture, &SrcR, &DestR);
  }
  }


  SDL_RenderPresent(s->renderer);
  s->nFrames++;
}

void * liveview_t(void * conf)
{
  liveXLview * xlconf = (liveXLview * ) conf;
  liveview(xlconf->X, xlconf->L, xlconf->N, &xlconf->quit);
  return NULL;
}

int liveview(double * X, uint8_t * L, size_t N, int * quit)
{

  printf("%d\n", __LINE__); fflush(stdout);

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

  sprintf(s->title, "live X view");

  // Initialize window
  gInit(s);
  getError(s);

  // Initialize beads

  bead * beads = malloc(26*sizeof(bead));
  for(int bb = 0; bb<26; bb++)
  {
    bead_init(s, beads+bb, bb);
  }

    getError(s);

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
    bead_free(s, &beads[bb]);
  }



  free(s->XZ);
  free(s->title);
  free(s);


  return 0;
}
