#include "forward.h"
#include <stdlib.h>

void forward(HMM *hmm, const unsigned int *Y, const unsigned int T, double * scalingFactor, double * alpha){
    
    unsigned int i;
    unsigned int j;
    unsigned int l;
    
    // Fill out the first row of alpha, using the init. probs.
    for(i = 0; i < hmm->hiddenStates; i++){
        alpha[i] = hmm->initProbs[i]*hmm->emissionProbs[i*hmm->observations+Y[0]];
        scalingFactor[0] += alpha[i];
    }
    
    // Scale the first row
    for(j = 0; j < hmm->hiddenStates; j++){
        alpha[j] = alpha[j]/scalingFactor[0];
    }
    
    // Fill out the rest of the rows in alpha
    for(i = 1; i < T; i++){
        for(j = 0; j < hmm->hiddenStates; j++){
            double emissionProb = hmm->emissionProbs[j*hmm->observations+Y[i]]; // Probability of emitting Y_i from state j
            
            double pastTransProb = 0.0; 
            for(l = 0; l < hmm->hiddenStates; l++){  // Emit the prefix, and transition K possible ways to the j'th state. Summed.
                pastTransProb += hmm->transitionProbs[l*hmm->hiddenStates+j]*alpha[(i-1)*hmm->hiddenStates+l];
            }
            
            alpha[i*hmm->hiddenStates+j] = emissionProb*pastTransProb; // Fill the i,j entry in alpha
            
            scalingFactor[i] += alpha[i*hmm->hiddenStates+j]; 
        }
        
        // Scale the i'th column of alpha
        for(j = 0; j < hmm->hiddenStates; j++){
            alpha[i*hmm->hiddenStates+j] = alpha[i*hmm->hiddenStates+j]/scalingFactor[i];
        }
    }
}
