#include "backward.h"
#include <stdlib.h>

void backward(HMM *hmm, const unsigned int *Y, const unsigned int T, double * scalingFactor, double * beta){
    
    unsigned int i;
    unsigned int j;
    
    // 2D beta matrix
    //
    // [state][time]
    //
    for(i = 0; i < hmm->hiddenStates; i++){
        beta[T*hmm->hiddenStates-1-i] = 1.0;
    }
    
    // Now for the step"BACKWARD" step
    for(i = T-1; i-- >0;){
        for(j = 0; j < hmm->hiddenStates; j++){
            for(int l = 0; l < hmm->hiddenStates; l++){
                double emissionProb = hmm->emissionProbs[l*hmm->observations+Y[i+1]];
                double transitionProb = hmm->transitionProbs[j*hmm->hiddenStates+l];
                double oldBeta = beta[(i+1)*hmm->hiddenStates+l];
                
                beta[i*hmm->hiddenStates+j] += transitionProb*emissionProb*oldBeta;
            }
            beta[i*hmm->hiddenStates+j] = beta[i*hmm->hiddenStates+j] / scalingFactor[i+1];
        }
    }
    
//    printf("Backward\n");
//    for(i = 0; i < T; i++){
//       for(j = 0; j < hmm->hiddenStates; j++){
//           printf("%f, ", beta[i*hmm->hiddenStates+j]);
//       }
//       printf("\n");
//    }
//    printf("\n");

    
}
