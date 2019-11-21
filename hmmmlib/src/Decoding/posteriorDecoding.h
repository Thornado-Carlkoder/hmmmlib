#ifndef posteriorDecoding_h
#define posteriorDecoding_h

#pragma once

#include <stdio.h>
#include "backward.h"
#include "hmm.h"

void posteriorDecoding(HMM * hmm, const unsigned int *Y, const int T, int * states);

#endif
