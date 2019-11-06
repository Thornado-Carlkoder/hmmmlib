#include "forward.h"
#include "hmm.h"
#include <stdlib.h>

void forward(HMM *hmm, const unsigned int *Y, const unsigned int T, double * scalingFactor, double * alpha){
    
    unsigned int i;
    unsigned int j;
    
    // 2D alpha matrix
    //
    // [state][time]
    //
    
    for(i = 0; i < hmm->hiddenStates; i++){
        alpha[i] = hmm->initProbs[i]*hmm->emissionProbs[i*hmm->observations+Y[0]];
        scalingFactor[0] += alpha[i];
    }
    
    // Scaling step
    for(j = 0; j < hmm->hiddenStates; j++){
        alpha[j] = alpha[j]/scalingFactor[0];
    }
    
    // Now the "recursive" step starts
    for(i = 1; i < T; i++){
        for(j = 0; j < hmm->hiddenStates; j++){
            double emissionProb = hmm->emissionProbs[j*hmm->observations+Y[i]];
            if(emissionProb > 0){
                double pastTransProb = 0.0;
                for(int l = 0; l < hmm->hiddenStates; l++){
                    pastTransProb += hmm->transitionProbs[l*hmm->hiddenStates+j]*alpha[(i-1)*hmm->hiddenStates+l];
                }
                alpha[i*hmm->hiddenStates+j] = emissionProb*pastTransProb;
            }
            scalingFactor[i] += alpha[i*hmm->hiddenStates+j];
        }
        // Scaling step
        for(j = 0; j < hmm->hiddenStates; j++){
            alpha[i*hmm->hiddenStates+j] = alpha[i*hmm->hiddenStates+j]/scalingFactor[i];
        }
    }
    
//    printf("Forward\n");
//    for(i = 0; i < T; i++){
//        for(j = 0; j < hmm->hiddenStates; j++){
//            printf("%f, ", alpha[i*hmm->hiddenStates+j]);
//        }
//        printf("\n");
//    }
//    printf("\n");
    
    //return alpha;
}
