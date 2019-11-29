#include <stdio.h>
#include <stdbool.h>
#pragma once

struct HMM {
    
    unsigned int hiddenStates;
    unsigned int observations;
    double * transitionProbs;
    double * emissionProbs;
    double * initProbs;
    
    void (*forward)(struct HMM *hmm, const unsigned int *Y, const unsigned int T, double * scalingFactor, double * alpha);
    void (*backward)(struct HMM *hmm, const unsigned int *Y, const unsigned int T, double * scalingFactor, double * beta);
    
};

typedef struct HMM HMM;

HMM * HMMConventionalsparse(const unsigned int hiddenStates, const unsigned int observations);

HMM * HMMConventional(const unsigned int hiddenStates, const unsigned int observations);

HMM * HMMBLAS(const unsigned int hiddenStates, const unsigned int observations);

HMM * HMMCsr(const unsigned int hiddenStates, const unsigned int observations);

HMM * HMMSBLAS(const unsigned int hiddenStates, const unsigned int observations);

bool validateHMM(const HMM *hmm);

void F(HMM *hmm, const unsigned int *Y, const unsigned int T, double * scalingFactor, double * alpha);

void B(HMM *hmm, const unsigned int *Y, const unsigned int T, double * scalingFactor, double * beta);

void printHMM(const HMM *hmm);

void HMMDeallocate(HMM * hmm);
