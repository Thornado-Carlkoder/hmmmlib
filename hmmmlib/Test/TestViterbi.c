#include "TestViterbi.h"


#include "viterbi.h"
#include "hmm.h"


#include <stdlib.h>
#include <assert.h>



extern bool testViterbi() {
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
        for(unsigned int j = 0; j < hmm->hiddenStates; j++){
            hmm->transitionProbs[i*hmm->hiddenStates+j] = transitionProbs[i][j];
        }
    }
    for(unsigned int i = 0; i < hmm->hiddenStates; i++){
        for(unsigned int j = 0; j < hmm->observations; j++){
            hmm->emissionProbs[i*hmm->observations+j] = emissionProbs[i][j];
        }
    }
    
    for(unsigned int i = 0; i < hmm->hiddenStates; i++){
        hmm->initProbs[i] = initProbs[i];
    }
    
    
    // obstest should be unsigned
    unsigned int obsTest1[99] = {3, 0, 2, 0, 2, 3, 3, 2, 3, 3, 2, 1, 3, 1, 0, 2, 3, 1, 0, 1, 0, 1, 0, 1, 3, 0, 1, 0, 2, 1, 0, 3, 0, 0, 0, 2, 3, 2, 0, 3, 0, 2, 3, 0, 2, 3, 0, 2, 2, 0, 1, 0, 2, 1, 0, 0, 3, 1, 2, 3, 1, 0, 3, 1, 0, 3, 1, 2, 0, 1, 0, 2, 0, 2, 1, 3, 0, 1, 2, 0, 3, 1, 0, 2, 3, 0, 2, 0, 0, 2, 0, 1, 0, 3, 1, 0, 3, 0, 0};
    //                 'T, A, G, A, G, T, T, G, T, T, G, C, T, C, A, G, T, C, A, C, A, C, A, C, T, A, C, A, G, C, A, T, A, A, A, G, T, G, A, T, A, G, T, A, G, T, A, G, G, A, C, A, G, C, A, A, T, C, G, T, C, A, T, C, A, T, C, G, A, C, A, G, A, G, C, T, A, C, G, A, T, C, A, G, T, A, G, A, A, G, A, C, A, T, C, A, T, A, A'
    
    
    unsigned int expectedTest1[99] = {3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1};

    unsigned int * obs = obsTest1;
    int nData = 99;
    
    // allocate matrix[data][states]
    int * table = calloc(nData*hmm->hiddenStates, sizeof(int)); // freed before return in viterbi()

    viterbi(hmm, obs, nData, table);
    
    
    for (int i = 0; i < nData; i++) {
        //printf("%u/%u ", viterbiResult[i], expectedTest1[i]);
        assert(table[i] == expectedTest1[i]);
    } 
    free(table);
    
    assert(valdidateHMM(hmm) == true);
    HMMDeallocate(hmm);
    
    return true;
}
