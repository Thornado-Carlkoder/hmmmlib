#include "hmm.h"
#pragma once

void backward_con_sparse(HMM *hmm, const unsigned int *Y, const unsigned int T, double * scalingFactor, double * beta);
