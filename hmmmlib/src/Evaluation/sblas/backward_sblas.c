#include "backward_sblas.h"
#include <stdlib.h>
#include <rsb.h>
//#include <Accelerate/Accelerate.h> // for mac os
#include <cblas.h> // for GNUlinux

void backward_sblas(HMM *hmm, const unsigned int *Y, const unsigned int T, double * scalingFactor, double * beta){
    
    unsigned int i;
    unsigned int j;
    
    //Creating the emission probs m*n into m, n*n matrix
    double ** new_emission_probs = calloc(hmm->observations, sizeof(double *));
    double * matrix = calloc(hmm->hiddenStates*hmm->hiddenStates, sizeof(double));
    
    for(i = 0; i < hmm->observations; i++){
        
        for(j = 0; j < hmm->hiddenStates; j++){
            matrix[j*hmm->hiddenStates+j] = hmm->emissionProbs[j*hmm->observations+i];
        }
        double * emission_probs = calloc(hmm->hiddenStates*hmm->hiddenStates, sizeof(double));
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, hmm->hiddenStates, hmm->hiddenStates, hmm->hiddenStates, 1.0, hmm->transitionProbs, hmm->hiddenStates, matrix, hmm->hiddenStates, 0.0, emission_probs, hmm->hiddenStates);
        new_emission_probs[i] = emission_probs;
    }
    free(matrix);
    
    const int bs = RSB_DEFAULT_BLOCKING;
    const int brA = bs, bcA = bs;
    const RSB_DEFAULT_TYPE one = 1;
    rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
    rsb_err_t errval = RSB_ERR_NO_ERROR;
    if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS)) != RSB_ERR_NO_ERROR)
    {
        printf("Error initializing the library!\n");
    }
    struct rsb_mtx_t ** mtx = malloc(hmm->observations*sizeof(struct rsb_mtx_t *));
    emission_rsb_mtx(hmm, new_emission_probs, mtx, &errval);

    for(i = 0; i < hmm->observations; i++){
        free(new_emission_probs[i]);
    }
    free(new_emission_probs);
    
    for(i = 0; i < hmm->hiddenStates; i++){
        beta[hmm->hiddenStates*T-1-i] = 1;
    }
    
    for(i = 1; i < T; i++){
        rsb_spmv(RSB_TRANSPOSITION_C, &one, mtx[Y[T-i]], beta+hmm->hiddenStates*T-i*hmm->hiddenStates, 1, &one, beta+hmm->hiddenStates*T-i*hmm->hiddenStates-hmm->hiddenStates, 1);
        cblas_dscal(hmm->hiddenStates, scalingFactor[T-i], beta+hmm->hiddenStates*T-i*hmm->hiddenStates-hmm->hiddenStates, 1);
    }
    
    for(i = 0; i < hmm->observations; i++){
        free(mtx[i]);
    }
    free(mtx);
}
