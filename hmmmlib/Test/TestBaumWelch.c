#include "TestBaumWelch.h"
#include "baumWelch.h"
#include "hmm.h"
#include <assert.h>
#include <math.h>

bool testBaumWelch() {
    HMM * hmm = HMMConventional(2, 2);
    
    double transitionProbs[2][2] = {
        {0.5, 0.5},
        {0.3, 0.7}
    };
    
    double emissionProbs[2][2] = {
        {0.3, 0.7},
        {0.8, 0.2}
    };
    
    double initProbs[2] = {0.2, 0.8};
    
    for(unsigned int i = 0; i < hmm->hiddenStates; i++){
        hmm->initProbs[i] = initProbs[i];
    }
    
    int i;
    int j;
    for(i = 0; i < hmm->hiddenStates; i++){
        for(j = 0; j < hmm->hiddenStates; j++){
            hmm->transitionProbs[i*hmm->hiddenStates+j] = transitionProbs[i][j];
        }
    }
    for(i = 0; i < hmm->hiddenStates; i++){
        for(j = 0; j < hmm->observations; j++){
            hmm->emissionProbs[i*hmm->observations+j] = emissionProbs[i][j];
        }
    }
    
    const unsigned int observation[10] = {0, 0, 0, 0, 0, 1, 1, 0, 0, 0};
    const unsigned int obsLenght = 10;
    
    baumWelch(hmm, observation, obsLenght, 1);
    
    // Checking the trans probs
    /*
     for(i = 0; i < hmm->hiddenStates*hmm->hiddenStates; i++){
     printf("%f, ", hmm->transitionProbs[i]);
     }
     */
    float epsilon = 0.000001;
    assert(fabs(hmm->transitionProbs[0]-0.43921478415604948) < epsilon);
    assert(fabs(hmm->transitionProbs[1*hmm->hiddenStates+1]-0.78554317773503979) < epsilon);
    
    // Checking the emission probs
    assert(fabs(hmm->emissionProbs[0]-0.46160107308583781) < epsilon);
    assert(fabs(hmm->emissionProbs[1*hmm->observations+1]-0.084984433203479412) < epsilon);
    
    //printHMM(hmm);
    
    assert(validateHMM(hmm) == true);
    HMMDeallocate(hmm);
    
    return true;

}
