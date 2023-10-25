#pragma once
/*
 * Conversion back and forth between HSV and RGB
 * See Alvy Ray Smith, Color Gamut Transform Pairs, SIGGRAPH '78.
 */

// Assumes that all values are in the range [0,1]

// Convert HSV to RGB. RGB will changed
void hsv2rgb(const double * restrict HSV, double * restrict RGB);

// Convert RGB to HSV. HSV will changed
void rgb2hsv(const double * restrict RGB, double * restrict HSV);
