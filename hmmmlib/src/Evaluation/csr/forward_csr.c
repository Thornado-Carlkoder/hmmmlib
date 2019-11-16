#include "forward_csr.h"
#include <stdlib.h>
#include <Accelerate/Accelerate.h> // for mac os
//#include <cblas.h> // for GNUlinux

int csr(HMM * hmm, double ** sparseMatrixs, int * ia, int * ja, double * a){
    
    unsigned int i;
    unsigned int j;
    unsigned int l;
    
    int nnz = 0;
    int nnz_count = 0;
    
    for(l = 0; l < hmm->observations; l++){
        for (i = 0; i < hmm->hiddenStates; i++) {
            for (j = 0; j < hmm->hiddenStates; j++) {
                if (sparseMatrixs[l][j*hmm->hiddenStates+i] != 0.0) {
                    a[nnz] = sparseMatrixs[l][j*hmm->hiddenStates+i];
                    if(l == 0){
                        ja[nnz] = j;
                        nnz_count++;
                    }
                    nnz++;
                }
            }
            if(l == 0){
                ia[i+1] = nnz;
            }
        }
    }
    return nnz_count;
}

void forward_csr(HMM *hmm, const unsigned int *Y, const unsigned int T, double * scalingFactor, double * alpha){

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
        
//        printf("------- Y[i] == %d -------\n", i);
//        for(int k = 0; k < hmm->hiddenStates; k++){
//            for(j = 0; j < hmm->hiddenStates; j++){
//                printf("%f, ", emission_probs[k*hmm->hiddenStates+j]);
//            }
//            printf("\n");
//        }
//        printf("\n\n\n");
        
    }
    free(matrix);
    
    // Doing the matrix multiplication and then scalingFactor
    scalingFactor[0] = cblas_dasum(hmm->hiddenStates, alpha, 1);
    cblas_dscal(hmm->hiddenStates, (1.0/scalingFactor[0]), alpha, 1);
    
    int * ia = calloc(hmm->hiddenStates+1, sizeof(int));
    int * ja = calloc(hmm->hiddenStates*hmm->hiddenStates, sizeof(int));
    double * a = calloc(hmm->observations*hmm->hiddenStates*hmm->hiddenStates, sizeof(double));
    
    int znn = csr(hmm, new_emission_probs, ia, ja, a);
    
    //printf("%d", znn);
    
    for(i = 0; i < hmm->observations; i++){
        free(new_emission_probs[i]);
    }
    free(new_emission_probs);
    
//    printf("\nAI: \n");
//    for(i = 0; i < hmm->hiddenStates+1; i++){
//        printf("%d, ", ia[i]);
//    }
//    printf("\nJA: \n");
//    for(i = 0; i < hmm->hiddenStates*hmm->hiddenStates; i++){
//        if(a[i] == 0.0) break;
//        printf("%d, ", ja[i]);
//    }
//
//    printf("\n\n");
    
    for(i = 1; i<T; i++){
        
        for (int row=0; row < hmm->hiddenStates; row++) {
            
            double sum = 0.0;
            
            for(int idx=ia[row]; idx<ia[row+1]; idx++) {
            
                sum += a[Y[i]*znn+idx] * alpha[(i-1)*hmm->hiddenStates+ja[idx]];
            
            }
            alpha[i*hmm->hiddenStates+row] = sum;
        }
        
        scalingFactor[i] = 1.0/cblas_dasum(hmm->hiddenStates, alpha+hmm->hiddenStates*i, 1);
        cblas_dscal(hmm->hiddenStates, scalingFactor[i], alpha+hmm->hiddenStates*i, 1);
    }
    
    free(ia);
    free(ja);
    free(a);

}
