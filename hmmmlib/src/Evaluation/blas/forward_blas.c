#include "forward_blas.h"
#include <stdlib.h>
//#include <Accelerate/Accelerate.h> // for mac os
#include <cblas.h> // for GNUlinux

void forward_blas(HMM *hmm, const unsigned int *Y, const unsigned int T, double * scalingFactor, double * alpha){
    
    unsigned int i;
    unsigned int j;
    
    //Creating the emission probs m*n into m, n*n matrix
    double ** new_emission_probs = calloc(hmm->observations, sizeof(double *));
    double * matrix = calloc(hmm->hiddenStates*hmm->hiddenStates, sizeof(double));
    
    for(i = 0; i < hmm->observations; i++){
        
        for(j = 0; j < hmm->hiddenStates; j++){
            matrix[j*hmm->hiddenStates+j] = hmm->emissionProbs[j*hmm->observations+i];
        }
        
        if(i == Y[0]){
            cblas_dsymv(CblasRowMajor, 121, hmm->hiddenStates, 1.0, matrix, hmm->hiddenStates, hmm->initProbs, 1, 1, alpha, 1);
        }
        
        double * emission_probs = calloc(hmm->hiddenStates*hmm->hiddenStates, sizeof(double));
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, hmm->hiddenStates, hmm->hiddenStates, hmm->hiddenStates, 1.0, hmm->transitionProbs, hmm->hiddenStates, matrix, hmm->hiddenStates, 0.0, emission_probs, hmm->hiddenStates);
        new_emission_probs[i] = emission_probs;
    }
    free(matrix);
    

    // Doing the matrix multiplication and then scalingFactor
    scalingFactor[0] = cblas_dasum(hmm->hiddenStates, alpha, 1);
    cblas_dscal(hmm->hiddenStates, (1.0/scalingFactor[0]), alpha, 1);

    for(i = 1; i<T; i++){

        cblas_dgemv(CblasRowMajor, CblasTrans, hmm->hiddenStates, hmm->hiddenStates, 1.0, new_emission_probs[Y[i]], hmm->hiddenStates, alpha+hmm->hiddenStates*(i-1), 1, 0, alpha+hmm->hiddenStates*i, 1);
        scalingFactor[i] = 1.0/cblas_dasum(hmm->hiddenStates, alpha+hmm->hiddenStates*i, 1);
        cblas_dscal(hmm->hiddenStates, scalingFactor[i], alpha+hmm->hiddenStates*i, 1);
        
    }
    for(i = 0; i < hmm->observations; i++){
        free(new_emission_probs[i]);
    }
    free(new_emission_probs);
}
