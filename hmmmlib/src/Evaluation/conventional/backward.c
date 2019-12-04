#include "backward.h"
#include <stdlib.h>

void backward(HMM *hmm, const unsigned int *Y, const unsigned int T, double * scalingFactor, double * beta){
    
    unsigned int i;
    unsigned int j;
    
    for(i = 0; i < hmm->hiddenStates; i++){
        beta[T*hmm->hiddenStates-1-i] = 1.0;
    }
    
    for(i = T-1; i-- >0;){
        for(j = 0; j < hmm->hiddenStates; j++){
            for(int l = 0; l < hmm->hiddenStates; l++){
                beta[i*hmm->hiddenStates+j] += hmm->transitionProbs[j*hmm->hiddenStates+l]*hmm->emissionProbs[l*hmm->observations+Y[i+1]]*beta[(i+1)*hmm->hiddenStates+l];
            }
            beta[i*hmm->hiddenStates+j] = beta[i*hmm->hiddenStates+j] / scalingFactor[i+1];
        }
    }
}
