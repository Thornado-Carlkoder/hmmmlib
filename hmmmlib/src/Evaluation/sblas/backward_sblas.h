#include "hmm.h"
#pragma once
void backward_sblas(HMM *hmm, const unsigned int *Y, const unsigned int T, double * scalingFactor, double * beta);

