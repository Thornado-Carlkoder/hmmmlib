#include "TestPosteriorDecoding.h"
#include "posteriorDecoding.h"
#include "hmm.h"
#include <stdlib.h>
#include <assert.h>

// skal linjen hvor "TestPosteriorDecoding.h" inkluderes ikke slettes?


extern bool testPosteriorDecoding() {
    HMM * hmm = HMMConventional(7, 4);
    
    double initProbs[7] = {0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00};
    
    double transitionProbs[7][7] ={
        {0.00, 0.00, 0.90, 0.10, 0.00, 0.00, 0.00},
        {1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00},
        {0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00},
        {0.00, 0.00, 0.05, 0.90, 0.05, 0.00, 0.00},
        {0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00},
        {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00},
        {0.00, 0.00, 0.00, 0.10, 0.90, 0.00, 0.00},
    };
    
    double emissionProbs[7][4] = {
        {0.30, 0.25, 0.25, 0.20},
        {0.20, 0.35, 0.15, 0.30},
        {0.40, 0.15, 0.20, 0.25},
        {0.25, 0.25, 0.25, 0.25},
        {0.20, 0.40, 0.30, 0.10},
        {0.30, 0.20, 0.30, 0.20},
        {0.15, 0.30, 0.20, 0.35}
    };


    for(unsigned int i = 0; i < hmm->hiddenStates; i++){
        hmm->initProbs[i] = initProbs[i];
    }
    for(unsigned int i = 0; i < hmm->hiddenStates; i++){
        for(unsigned int j = 0; j < hmm->hiddenStates; j++){
            hmm->transitionProbs[i*hmm->hiddenStates+j] = transitionProbs[i][j];
        }
    }
    
    for(unsigned int i = 0; i < hmm->hiddenStates; i++){
        for(unsigned int j = 0; j < hmm->observations; j++){
            hmm->emissionProbs[i*hmm->observations+j] = emissionProbs[i][j];
        }
    }
    

    unsigned int data[105] =     {0, 1, 3, 3, 2, 1, 2, 3, 0, 0, 1, 2, 1, 0, 2, 3, 0, 1, 3, 1, 2, 0, 3, 1, 2, 0, 3, 1, 2, 3, 0, 1, 2, 1, 0, 3, 2, 1, 0, 0, 1, 0, 3, 2, 0, 3, 0, 0, 0, 3, 1, 1, 0, 0, 3, 1, 0, 3, 1, 0, 3, 0, 1, 0, 3, 2, 0, 1, 0, 3, 1, 0, 2, 1, 3, 0, 2, 1, 0, 3, 2, 0, 1, 2, 0, 3, 1, 0, 2, 3, 0, 1, 2, 0, 1, 3, 0, 1, 2, 0, 3, 1, 2, 3, 0};
    unsigned int expected[105] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0};    
    unsigned int data_len = 105;

    unsigned int * states = calloc(data_len, sizeof(int));
    
    posteriorDecoding(hmm, data, data_len, states);

    // debug print output from algorithm
    for (unsigned int i = 0; i < data_len; i++) {
        //printf("%u(%u) ", result[i], expected[i]);
        assert(states[i] == expected[i]);
    }
    
    free(states);
    assert(validateHMM(hmm) == true);
    HMMDeallocate(hmm);

    return true;
}
