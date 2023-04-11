#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

static double max(double a, double b)
{
  if(a>b)
  {
    return a;
  }
  return b;
}


static double min(const double a, const double b)
{
  if(a<b)
  {
    return a;
  }
  return b;
}

void hsv2rgb(const double * restrict HSV, double * restrict RGB)
{
    // 1
  double H = HSV[0]*6;
  double S = HSV[1];
  double V = HSV[2];
// 2

  if(H==6)
    H = 0;

  int I = floor(H);

  double F = H-I;
  // 3
  double M = V*(1-S);
  double N = V*(1-S*F);
  double K = V*(1-S*(1-F));
  // 4
  
//  printf("I: %d\n", I);

  switch(I)
  {
    case 0:
      RGB[0] = V; RGB[1] = K; RGB[2] = M;
      break;
    case 1:
      RGB[0] = N; RGB[1] = V; RGB[2] = M;
      break;
    case 2:
      RGB[0] = M; RGB[1] = V; RGB[2] = K;
      break;
    case 3:
      RGB[0] = M; RGB[1] = N; RGB[2] = V;
      break;
    case 4:
      RGB[0] = K; RGB[1] = M; RGB[2] = V;
      break;
    case 5:
      RGB[0] = V; RGB[1] = M; RGB[2] = N;
      break;
  }

return;
}

void rgb2hsv(const double * restrict RGB, double * restrict HSV)
{

  double H = 0; double S = 0; double V = 0;

  double R = RGB[0];
  double G = RGB[1];
  double B = RGB[2];

  // 1
  V = max(R, max(G, B));
  // 2
  double X = min(R, min(G, B));
  // 3
  S = (V-X)/V;
  if(S == 0 || !isfinite(S))
  {      
    HSV[0] = 0; HSV[1] = 0; HSV[2] = V;
    return;
  }
  // 4
  double r = (V-R)/(V-X);
  double g = (V-G)/(V-X);
  double b = (V-B)/(V-X);
  // 5
 // printf("R %f G %f B %f r %f g %f b %f\n", R, G, B, r, g, b);
  if(R == V)
  {
    if(G == X)
    {
      H = 5+b;
    } else {
      H = 1-g;
    }
  }

  if(G == V)
  {
    if(B == X)
    {
      H = 1+r;
    } else {
      H = 3-b;
    }
  }

  if(B == V) 
  {
    if(R == X)
    {
      H = 3+g;
    } else {
      H = 5-r;
    }
  }

  // 6
  H = H/6;

  //
  HSV[0] = H; HSV[1] = S; HSV[2] = V;
  return;
}


