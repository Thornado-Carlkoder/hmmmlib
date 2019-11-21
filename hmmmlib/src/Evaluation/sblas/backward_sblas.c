#include "backward_sblas.h"
#include <stdlib.h>
#include <rsb.h>
//#include <Accelerate/Accelerate.h> // for mac os
#include <cblas.h> // for GNUlinux

void b_emission_rsb_mtx(HMM * hmm, double ** sparseMatrixs, struct rsb_mtx_t ** rsb_mtx, rsb_err_t * errval){
    
    rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
    const int brA = RSB_DEFAULT_BLOCKING;
    const int bcA = RSB_DEFAULT_BLOCKING;
    
    unsigned int i;
    unsigned int j;
    
    int * ia = calloc(hmm->hiddenStates*hmm->hiddenStates, sizeof(int));
    int * ja = calloc(hmm->hiddenStates*hmm->hiddenStates, sizeof(int));
    
    int nnz = 0;
    
    for (j = 0; j < hmm->hiddenStates; j++) {
        for (i = 0; i < hmm->hiddenStates; i++) {
            if (sparseMatrixs[0][j*hmm->hiddenStates+i] != 0.0) {
                ja[nnz] = i;
                ia[nnz] = j;
                nnz++;
            }
        }
    }
    
    int * final_ja = malloc(nnz*sizeof(int));
    int * final_ia = malloc(nnz*sizeof(int));
    double * final_a = malloc(nnz*sizeof(double));
    
    for(i = 0; i < nnz; i++){
        final_ja[i] = ja[i];
        final_ia[i] = ia[i];
    }
    
    free(ia);
    free(ja);

    for(i = 0; i < hmm->observations; i++){
        for(j = 0; j < nnz; j++){
            final_a[j] = sparseMatrixs[i][final_ia[j]*hmm->hiddenStates+final_ja[j]];
            //printf("%f, ", final_a[j]);
        }
        rsb_mtx[i] = rsb_mtx_alloc_from_coo_const(final_a, final_ia, final_ja, nnz, typecode, hmm->hiddenStates, hmm->hiddenStates, brA, bcA, RSB_FLAG_NOFLAGS | RSB_FLAG_DUPLICATES_SUM, errval);
    }
    
    free(final_ja);
    free(final_ia);
    free(final_a);
}

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
    b_emission_rsb_mtx(hmm, new_emission_probs, mtx, &errval);

    for(i = 0; i < hmm->observations; i++){
        free(new_emission_probs[i]);
    }
    free(new_emission_probs);
    
    for(i = 0; i < hmm->hiddenStates; i++){
        beta[hmm->hiddenStates*T-1-i] = 1;
    }
    
    for(i = 1; i < T; i++){
        rsb_spmv(RSB_TRANSPOSITION_N, &one, mtx[Y[T-i]], beta+hmm->hiddenStates*T-i*hmm->hiddenStates, 1, &one, beta+hmm->hiddenStates*T-i*hmm->hiddenStates-hmm->hiddenStates, 1);
        cblas_dscal(hmm->hiddenStates, scalingFactor[T-i], beta+hmm->hiddenStates*T-i*hmm->hiddenStates-hmm->hiddenStates, 1);
    }

    for(i = 0; i < hmm->observations; i++){
        free(mtx[i]);
    }
    free(mtx);
}
