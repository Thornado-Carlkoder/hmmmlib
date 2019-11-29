#include "hmm.h"
#include <stdlib.h>
#include <time.h>

void runningTimeTest(){
    HMM * hmm2 = HMMCsr(7, 4);
    
    int i;
    int j;
    
    double transitionProbs2[7][7] = {
     {0.0 , 0.0 , 0.9 , 0.1 , 0.0 , 0.0 , 0.0},
     {1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0},
     {0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0},
     {0.0 , 0.0 , 0.05 , 0.9 , 0.05 , 0.0 , 0.0},
     {0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0},
     {0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0},
     {0.0 , 0.0 , 0.0 , 0.1 , 0.9 , 0.0 , 0.0}
    };

    double emissionProbs2[7][4] = {
     {0.3 , 0.25 , 0.25 , 0.2},
     {0.2 , 0.35 , 0.15 , 0.3},
     {0.4 , 0.15 , 0.2 , 0.25},
     {0.25 , 0.25 , 0.25 , 0.25},
     {0.2 , 0.4 , 0.3 , 0.1},
     {0.3 , 0.2 , 0.3 , 0.2},
     {0.15 , 0.3 , 0.2 , 0.35}
    };

    double initProbs2[7] = {0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0};

    for(i = 0; i < hmm2->hiddenStates; i++){
      hmm2->initProbs[i] = initProbs2[i];
    }
    for(i = 0; i < hmm2->hiddenStates; i++){
     for(j = 0; j < hmm2->hiddenStates; j++){
         hmm2->transitionProbs[i*hmm2->hiddenStates+j] = transitionProbs2[i][j];
     }
    }
    for(i = 0; i < hmm2->hiddenStates; i++){
     for(j = 0; j < hmm2->observations; j++){
         hmm2->emissionProbs[i*hmm2->observations+j] = emissionProbs2[i][j];
     }
    }
    
    const unsigned int obsLenght2 = 1000000;
    unsigned int observation2[obsLenght2];
    srand(time(NULL));
    for(i = 0; i < obsLenght2; i++){
        observation2[i] = rand() % 4;
    }
    clock_t t;
    double time_taken;
    double * scaleFactor2 = calloc(obsLenght2, sizeof(double));
    double * alpha2 = calloc(obsLenght2*hmm2->hiddenStates, sizeof(double));
    double * beta2 = calloc(obsLenght2*hmm2->hiddenStates, sizeof(double));
    for(i = 0; i < 10; i++){
        t = clock();
        F(hmm2, observation2, obsLenght2, scaleFactor2, alpha2);
        t = clock() - t;
        time_taken = ((double)t)/CLOCKS_PER_SEC;
        printf("Time Forward: 7x7 %dk obs: %f sec\n", obsLenght2/1000, time_taken);
    }
    for(i = 0; i < 10; i++){
        t = clock();
        B(hmm2, observation2, obsLenght2, scaleFactor2, beta2);
        t = clock() - t;
        time_taken = ((double)t)/CLOCKS_PER_SEC;
        printf("Time Backward: 7x7 %dk obs: %f sec\n", obsLenght2/1000, time_taken);
    }
    free(scaleFactor2);
    free(alpha2);
    free(beta2);
    HMMDeallocate(hmm2);
    printf("Running time blas!!!");
  
}
