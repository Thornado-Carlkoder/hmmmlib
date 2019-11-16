#include "forward_sblas.h"
#include <stdlib.h>
#include <rsb.h>
//#include <Accelerate/Accelerate.h> // for mac os
#include <cblas.h> // for GNUlinux

int hulla_csr(HMM * hmm, double ** sparseMatrixs, struct rsb_mtx_t ** rsb_mtx, rsb_err_t * errval){
    
    rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
    const int bs = RSB_DEFAULT_BLOCKING;
    const int brA = bs, bcA = bs;
    
    unsigned int i;
    unsigned int j;
    
    int * ia = calloc(hmm->hiddenStates*hmm->hiddenStates, sizeof(int));
    int * ja = calloc(hmm->hiddenStates*hmm->hiddenStates, sizeof(int));
    
    //double * a = calloc(hmm->observations*hmm->hiddenStates*hmm->hiddenStates, sizeof(double));
    //free(a)
    
    int nnz = 0;
    
    for (j = 0; j < hmm->hiddenStates; j++) {
        for (i = 0; i < hmm->hiddenStates; i++) {
            if (sparseMatrixs[0][j*hmm->hiddenStates+i] != 0.0) {
                //a[nnz] = sparseMatrixs[0][j*hmm->hiddenStates+i];
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
        }
        rsb_mtx[i] = rsb_mtx_alloc_from_coo_const(final_a, final_ia, final_ja, nnz, typecode, hmm->hiddenStates, hmm->hiddenStates, brA, bcA, RSB_FLAG_NOFLAGS | RSB_FLAG_DUPLICATES_SUM, errval);
    }
    
//    printf("\n\n------------------------\n");
//    for(i = 0; i<hmm->observations; i++){
//        for(j=0; j < nnz; j++){
//            printf("%f, ",a[i*nnz+j]);
//        }
//        printf("\n");
//    }
//    printf("\n------------------------\n\n");
    
    return nnz;
}

void forward_sblas(HMM *hmm, const unsigned int *Y, const unsigned int T, double * scalingFactor, double * alpha){
    
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
    
    const int bs = RSB_DEFAULT_BLOCKING;
    const int brA = bs, bcA = bs;
    const RSB_DEFAULT_TYPE one = 1;
    rsb_type_t typecode = RSB_NUMERICAL_TYPE_DEFAULT;
    rsb_err_t errval = RSB_ERR_NO_ERROR;
    struct rsb_mtx_t *mtxAp = NULL; /* matrix structure pointer */
    rsb_coo_idx_t IA[] = {0,1,1,2};
    /* nonzero column indices coordinates: */
    rsb_coo_idx_t JA[] = {0,1,2,2};
    RSB_DEFAULT_TYPE VA[] = {11.5,10.5,22.5,32.5};/* values of nonzeroes */

    printf("Hello, RSB!\n");
    printf("Initializing the library...\n");
    if((errval = rsb_lib_init(RSB_NULL_INIT_OPTIONS)) != RSB_ERR_NO_ERROR)
    {
        printf("Error initializing the library!\n");
    }
    struct rsb_mtx_t ** mtx = malloc(hmm->observations*sizeof(struct rsb_mtx_t *));
    int znn = hulla_csr(hmm, new_emission_probs, mtx, &errval);
    //int znn = 11;
    for(i = 0; i < hmm->observations; i++){
        free(new_emission_probs[i]);
    }
    free(new_emission_probs);
    
    mtxAp = rsb_mtx_alloc_from_coo_const(VA, IA, JA, znn, typecode, hmm->hiddenStates, hmm->hiddenStates, brA, bcA, RSB_FLAG_NOFLAGS | RSB_FLAG_DUPLICATES_SUM, &errval);

    for(i = 1; i<T; i++){
        rsb_spmv(RSB_TRANSPOSITION_C, &one, mtx[Y[i]], alpha+hmm->hiddenStates*(i-1), 1, &one, alpha+hmm->hiddenStates*i, 1);
        scalingFactor[i] = 1.0/cblas_dasum(hmm->hiddenStates, alpha+hmm->hiddenStates*i, 1);
        cblas_dscal(hmm->hiddenStates, scalingFactor[i], alpha+hmm->hiddenStates*i, 1);
    }
    
    for(i = 0; i < hmm->observations; i++){
        free(mtx[i]);
    }
    free(mtx);
}
