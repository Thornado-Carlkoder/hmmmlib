#include <stdio.h>
#include "hmm.h"

#pragma once

void forward_con_sparse(HMM *hmm, const unsigned int *Y, const unsigned int T, double * scalingFactor, double * alpha);
