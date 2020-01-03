#include "forward_blas.h"
#include <stdlib.h>
//#include <Accelerate/Accelerate.h> // for mac os
#include <cblas.h> // for GNUlinux

void forward_blas(HMM *hmm, const unsigned int *Y, const unsigned int T, double * scalingFactor, double * alpha){
    
    unsigned int i;
    unsigned int j;
    double ** new_emission_probs = calloc(hmm->observations, sizeof(double *));
    double * matrix = calloc(hmm->hiddenStates*hmm->hiddenStates, sizeof(double)); // Zero matrix
    
    for(i = 0; i < hmm->observations; i++){ // For each symbol in the alphabet
        
        
        // A diagonal matrix 'matrix[j]' is created for each symbol in the alphabet
        // It contains the probability of emitting the i'th symbol from state j
        for(j = 0; j < hmm->hiddenStates; j++){ // (For each hidden state in the state space.)
            matrix[j*hmm->hiddenStates+j] = hmm->emissionProbs[j*hmm->observations+i]; 
        }
        
        
        // If the i'th symbol is equal to the first in the data:
        // Use BLAS to compute the following:
        // alpha :=  matrix * hmm->initProbs + alpha
        if(i == Y[0]){ 
            cblas_dsymv(CblasRowMajor, 121, hmm->hiddenStates,   1.0, matrix, hmm->hiddenStates, hmm->initProbs,    1,    1, alpha,    1);
            //          Layout        uplo,                 n, alpha,     *a,               lda,              x, incx, beta,     y, incy
        }
        
        double * emission_probs = calloc(hmm->hiddenStates*hmm->hiddenStates, sizeof(double));
        
        // Use BLAS to calculate the following:
        // emission_probs = hmm->transitionProbs * matrix      
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, hmm->hiddenStates, hmm->hiddenStates, hmm->hiddenStates,     1.0, hmm->transitionProbs, hmm->hiddenStates, matrix, hmm->hiddenStates,  0.0, emission_probs, hmm->hiddenStates);
        //                  Order,       TransA,       TransB,                 M,                 N,                 K,   alpha,                    A,               lda,      B                ldb, beta,              C,               ldc
        
        // Save this emission_probs in a list.
        new_emission_probs[i] = emission_probs;
    }
    free(matrix);
    
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
