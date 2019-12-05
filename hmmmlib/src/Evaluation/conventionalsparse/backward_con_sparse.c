#include "backward_con_sparse.h"
#include <stdlib.h>

void backward_con_sparse(HMM *hmm, const unsigned int *Y, const unsigned int T, double * scalingFactor, double * beta){
    
    unsigned int i;
    unsigned int j;
    unsigned int l;
    
    
    for(i = 0; i < hmm->hiddenStates; i++){
        beta[T*hmm->hiddenStates-1-i] = 1.0;
    }
    
    for(i = T-1; i-- >0;){
        for(j = 0; j < hmm->hiddenStates; j++){
            for(l = 0; l < hmm->hiddenStates; l++){
                double transitionProb = hmm->transitionProbs[j*hmm->hiddenStates+l];
                if(transitionProb > 0){
                    double emissionProb = hmm->emissionProbs[l*hmm->observations+Y[i+1]];
                    double oldBeta = beta[(i+1)*hmm->hiddenStates+l];
                    beta[i*hmm->hiddenStates+j] += transitionProb*emissionProb*oldBeta;
                }
            
            }
            beta[i*hmm->hiddenStates+j] = beta[i*hmm->hiddenStates+j] / scalingFactor[i+1];
        }
    }
}
