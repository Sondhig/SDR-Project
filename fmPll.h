/*
Comp Eng 3DY4 (Computer Systems Integration Project)
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_fmPll_H
#define DY4_fmPll_H

// add headers as needed
#include <cmath>
#include <math.h>

void fmPll(std::vector<float> &pllIn, int freq, float Fs,
           std::vector<float> &state, float ncoScale, float phaseAdjust,
           float normBandwidth, std::vector<float> &ncoOut);
#endif
