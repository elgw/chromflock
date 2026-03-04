#include "liveview.h"

typedef uint32_t u32;
typedef int64_t i64;
typedef uint8_t u8;

typedef struct{
    double X[3];
    u32 label;
} bead;

#define BEAD_WIDTH 51
#define BEAD_HEIGHT 51

typedef struct {
    size_t N; // number of beads

    const uint8_t * L; // labels
    const double * X; // pointer to "live" data

    bead * beads;

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
    const elli * E;
    double * EA;
    int perspective; // if 0: orthographic projection
    int show_lines;
} scene;

typedef struct {
    u8 r;
    u8 g;
    u8 b;
    u8 a;
} pixel;

typedef struct {
    SDL_Texture * texture;
} bead_graphics;


float mousex = 0;
float mousey = 0;
int mousedown = 0;

// T = transpose(A), T: 3x3
static void
mattrans(double * T, double * A)
{
    T[0] = A[0]; T[3] = A[1]; T[6] = A[2];
    T[1] = A[3]; T[4] = A[4]; T[7] = A[5];
    T[2] = A[6]; T[5] = A[7]; T[8] = A[8];
}

// In: X, Y : 3x3 matrices
// Out: Z = X*Y,
static void
matmul(double * Z, const double * X, const double * Y)
{
    double T[9] = {0};
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
    for(int kk = 0; kk < 9; kk++)
    {
        Z[kk] = T[kk];
    }
    return;
}

// In: R: 3x3, X: 3x1
// Out: X=R*X: 3x1
static void
rot_point(const double * restrict R, double * restrict X)
{
    double x = R[0] * X[0] + R[3]*X[1] + R[6]*X[2];
    double y = R[1] * X[0] + R[4]*X[1] + R[7]*X[2];
    double z = R[2] * X[0] + R[5]*X[1] + R[8]*X[2];
    X[0] = x;
    X[1] = y;
    X[2] = z;
}

// Rotate around x - axis
// SR = SR*Rx(theta)
static void
rot_x(double * SR, double theta)
{
    double R[9];
    R[0] = 1; R[3] = 0;           R[6] = 0;
    R[1] = 0; R[4] = cos(theta);  R[7] = sin(theta);
    R[2] = 0; R[5] = -sin(theta); R[8] = cos(theta);
    matmul(SR, R, SR);
}

// Rotate around y - axis
static void
rot_y(double * SR, const double theta)
{
    double R[9];
    R[0] = cos(theta);  R[3] = 0; R[6] = sin(theta);
    R[1] = 0;           R[4] = 1; R[7] = 0;
    R[2] = -sin(theta); R[5] = 0; R[8] = cos(theta);
    matmul(SR, R, SR);
}


// Rotate around z - axis
static void
rot_z(double * SR, const double theta)
{
    double R[9];
    R[0] =  cos(theta); R[3] = sin(theta); R[6] = 0;
    R[1] = -sin(theta); R[4] = cos(theta); R[7] = 0;
    R[2] = 0;           R[5] = 0;          R[8] = 1;
    matmul(SR, R, SR);
}


static double cmap[] = {255,255,255,
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

static void
drawBead(pixel* pixels,
         int width, int height, const int label)
{

    if(width != height)
    {
        printf("Width not equal to height\n");
        exit(0);
    }

    int color_id = label;
    color_id = color_id % 32;
    if(color_id > 24)
    {
        color_id = 0;
    }

    double base_hsv[4] = {0};
    double base_rgb[4] = {
        cmap[3*color_id+0]/255.0,
        cmap[3*color_id+1]/255.0,
        cmap[3*color_id+2]/255.0,
        1.0};

    rgb2hsv(base_rgb, base_hsv);

    int br = round((width-1)/2); // bead radius

    for(int xx = -br; xx<=br; xx++)
    {
        for(int yy = -br; yy<=br; yy++)
        {
            double rad = sqrt(pow(xx,2) + pow(yy, 2));
            i64 idx = (br+xx)+width*(br+yy); // pixel idx
            double pixel_hsv[4] = {0};
            double pixel_rgb[4] = {0};
            if(rad <= br)
            {
                memcpy(pixel_hsv, base_hsv, 4*sizeof(double));
                pixel_hsv[2] *= sqrt((br-rad)/br); // change "value"
                hsv2rgb(pixel_hsv, pixel_rgb);
                pixel_rgb[3] = 1.0;
                if(label >= 32)
                {
                    if(rad < br/5.0 )
                    {
                        pixel_rgb[0] = 1.0 - pixel_rgb[0];
                        pixel_rgb[1] = 1.0 - pixel_rgb[1];
                        pixel_rgb[2] = 1.0 - pixel_rgb[2];
                    }
                }
            } else {
                pixel_rgb[3] = 0;
            }


            if(0)
            {
                printf("%d, %d : %f %f %f %f\n",
                       xx, yy,
                       pixel_rgb[0],
                       pixel_rgb[1],
                       pixel_rgb[2],
                       pixel_rgb[3]);
            }
            pixels[idx].r = (u8) (pixel_rgb[0]*255.0);
            pixels[idx].g = (u8) (pixel_rgb[1]*255.0);
            pixels[idx].b = (u8) (pixel_rgb[2]*255.0);
            pixels[idx].a = (u8) (pixel_rgb[3]*255.0);
        }
    }

    return;
}

void
bead_init(scene * s, bead_graphics * b, int label)
{
    int width = BEAD_WIDTH;
    int height = BEAD_HEIGHT;
    int depth = 32;
    int pitch = 4*width;

    pixel * pixels = calloc(width*height, sizeof(pixel));
    assert(pixels != NULL);

    drawBead(pixels, width, height, label);

    u32 rmask = 0x000000ff;
    u32 gmask = 0x0000ff00;
    u32 bmask = 0x00ff0000;
    u32 amask = 0xff000000;

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

void bead_graphics_free(bead_graphics *b)
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

            if (evt.key.keysym.sym == SDLK_p)
            {
                s->perspective++;
                s->perspective = s->perspective % 2;
            }

            if (evt.key.keysym.sym == SDLK_l)
            {
                s->show_lines++;
                s->show_lines = s->show_lines % 2;
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
    bead * P = (bead * ) A;
    bead * Q = (bead * ) B;

    if(P->X[2] > Q->X[2])
        return 1;
    if(P->X[2] < Q->X[2])
        return -1;

    return 0;
}

// Make a copy of the input array XYZ
// and sort it by Z value after it is rotated
static void
copy_sort(scene * s)
{
    for(size_t kk = 0; kk < s->N; kk++)
    {
        s->beads[kk].X[0] = s->X[3*kk];
        s->beads[kk].X[1] = s->X[3*kk+1];
        s->beads[kk].X[2] = s->X[3*kk+2];
        s->beads[kk].label = s->L[kk];
        rot_point(s->Rot, s->beads[kk].X);
    }
    qsort(s->beads, s->N, sizeof(bead), zcmp);
    return;
}


static void
render(scene * s, const bead_graphics * beads)
{
    // Only needed if s->X is updated
    if(s->pause == 0)
    {
        copy_sort(s);
    }

    SDL_SetRenderDrawColor(s->renderer, 255.0, 255.0, 255.0, SDL_ALPHA_OPAQUE);
    SDL_RenderClear(s->renderer);

    // Bead radius
    const double bead_radius = s->r0;
    // In terms of screen pixels
    const int br = round((double) imin(s->window_w, s->window_h)*bead_radius);

    const int mid = imin(s->window_w / 2, s->window_h /2);

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

    // draw a circle for the domain
    // actually draws 3 circles to make it thicker
    // TODO: use SDL_RnederDrawLines
    for(int delta = 0; delta<3; delta++)
    {
        double x0 = -10;
        double y0 = -10;

        for(double theta = 0.01; theta<=6*M_PI; theta = theta+0.02)
        {
            double x = cos(theta); double y = sin(theta);
            double scale = // ||x^TAx||
                x*(ER[0]*x + ER[3]*y) + y*(ER[1]*x + ER[4]*y);

            x = x*(mid+delta)/sqrt(scale) + mid + woff;
            y = y*(mid+delta)/sqrt(scale) + mid + hoff;

            if(x0>-10)
                SDL_RenderDrawLine(s->renderer, round(x0), round(y0), round(x), round(y));
            y0 = y;
            x0 = x;
        }
    }

    if(s->show_lines == 1)
    {
        for(size_t kk = 0; kk+1 < s->N; kk++)
        {
            // quite silly without a z-buffer or similar
            // also silly with 1 px thick lines
            if(s->L[kk] != s->L[kk+1])
                continue;
            int label = s->L[kk] % 32;
            size_t i0 = 3*kk;
            size_t i1 = 3*(kk+1);
            const double * X = s->X;
            double X0[3] = {X[i0], X[i0+1], X[i0+2]};
            double X1[3] = {X[i1], X[i1+1], X[i1+2]};
            rot_point(s->Rot, X0);
            rot_point(s->Rot, X1);

            SDL_SetRenderDrawColor(s->renderer,
                                   (u8) cmap[3*label],
                                   (u8) cmap[3*label+1],
                                   (u8) cmap[3*label+2], 255);

            // project to screen
            X0[0] = mid + mid*X0[0] + woff;
            X1[0] = mid + mid*X1[0] + woff;
            X0[1] = mid + mid*X0[1] + hoff;
            X1[1] = mid + mid*X1[1] + hoff;
            //printf("%f, %f -- %f, %f\n", X0[0], X0[1], X1[0], X1[1]);
            SDL_RenderDrawLine(s->renderer,
                               round(X0[0]), round(X0[1]),
                               round(X1[0]), round(X1[1]));
        }
    }

    /* Draw beads */
    if(s->show_lines == 0)
    {
        for(size_t kk = 0; kk<s->N; kk++)
        {
            int label = (int) s->beads[kk].label;

            SDL_Rect SrcR;
            SDL_Rect DestR;

            SrcR.x = 0;
            SrcR.y = 0;
            SrcR.w = BEAD_WIDTH;
            SrcR.h = BEAD_HEIGHT;

            SDL_QueryTexture(beads[label].texture, NULL, NULL, &SrcR.w, &SrcR.h);

            DestR.w = br;
            DestR.h = br;
            double dest_x = s->beads[kk].X[0];
            double dest_y = s->beads[kk].X[1];

            if(s->perspective)
            {
                // Beads are sorted with the smallest z value first (-1)
                // and the largest z value last (1) which will be on the top
                const double p = 7.0;
                const double ps = p/(p - s->beads[kk].X[2]); // projective scaling
                dest_x  *= ps;
                dest_y  *= ps;
                DestR.w *= ps;
                DestR.h *= ps;
            }

            DestR.x = mid + mid*dest_x - DestR.w/2 + woff;
            DestR.y = mid + mid*dest_y - DestR.w/2 + hoff;

            SDL_RenderCopy(s->renderer, beads[label].texture, &SrcR, &DestR);
        }
    }

    SDL_RenderPresent(s->renderer);
    s->nFrames++;
}


int
liveview(const double * X,
         const uint8_t * L,
         const size_t N,
         volatile int * quit,
         const double r0,
         const elli * E)
{
    scene * s = calloc(1, sizeof(scene));
    assert(s != NULL);

    s->X = X;
    s->N = N;
    s->r0 = r0;
    s->beads = calloc(s->N, sizeof(bead));
    assert(s->XZ != NULL);
    s->X = X;
    s->L = L;
    s->done = 0;
    s->window_h = 512;
    s->window_w = 512;
    s->title = calloc(512, sizeof(char));
    assert(s->title != NULL);
    s->E = E;
    s->Rot = calloc(9, sizeof(double));
    assert(s->Rot != NULL);
    s->Rot[0] = 1; s->Rot[4] = 1; s->Rot[8] = 1;
    s->EA = calloc(9, sizeof(double));
    assert(s->EA != NULL);
    s->EA[0] = 1.0/pow(s->E->a, 2);
    s->EA[4] = 1.0/pow(s->E->b, 2);
    s->EA[8] = 1.0/pow(s->E->c, 2);
    sprintf(s->title, "live X view");

    // Initialize window
    gInit(s);

    // Initialize beads
    bead_graphics * tbeads = calloc(64, sizeof(bead_graphics));
    assert(beads != NULL);
    for(int bb = 0; bb<64; bb++)
    {
        bead_init(s, tbeads+bb, bb);
    }

    // main loop
    while(s->done == 0 && quit[0] == 0)
    {
        render(s, tbeads);
        getEvents(s);
        // usleep(1000000.0/60.0);
    }

    fprintf(stdout, "Rendered %zu times\n", s->nFrames);

    SDL_DestroyRenderer(s->renderer);
    SDL_DestroyWindow(s->window);

    for(int bb = 0; bb<64; bb++)
    {
        bead_graphics_free(&tbeads[bb]);
    }
    free(tbeads);
    free(s->beads);
    free(s->title);
    free(s);

    return EXIT_SUCCESS;
}
