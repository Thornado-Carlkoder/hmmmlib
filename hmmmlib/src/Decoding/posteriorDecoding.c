#include "posteriorDecoding.h"
#include <stdlib.h>
#include <math.h>

void posteriorDecoding(HMM * hmm, const unsigned int *Y, const unsigned int T, unsigned int * states){
    
    double * alpha = calloc(hmm->hiddenStates*T, sizeof(double));
    double * beta = calloc(hmm->hiddenStates*T, sizeof(double));
    double * scalingFactor = calloc(T, sizeof(double));
    
    F(hmm, Y, T, scalingFactor, alpha);
    B(hmm, Y, T, scalingFactor, beta);

    double * z_ = calloc(T, sizeof(double));
    
    //double posterior; // Skal posterior ikke være en liste?
    double * posterior = calloc(hmm->hiddenStates, sizeof(double));

    for (unsigned int i = 0; i < T; i++) {
        for (unsigned int k = 0; k < hmm->hiddenStates; k++) {
            //posterior[k] = alpha[i*T+k]*beta[i*T+k];
            posterior[k] = alpha[i*hmm->hiddenStates+k] * beta[i*hmm->hiddenStates+k];
        }

        // argmax posterior
        unsigned int maxIndex = 0; // initialization value
        double max = -INFINITY;
        for (unsigned int l = 0; l < hmm->hiddenStates; l++) {
            if (posterior[l] > max) {
                maxIndex = l;
                max = posterior[l];
                //printf("max is %f\n", max);
            }
        }
        states[i] = maxIndex;
        z_[i] = max;
        
    }
    
    free(scalingFactor);
    free(alpha);
    free(beta);
    free(posterior);
    free(z_);
}
