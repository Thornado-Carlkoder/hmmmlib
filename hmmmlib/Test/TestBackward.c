#include "TestBackward.h"
#include "hmm.h"

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

bool testBackward2x2(){
    
    HMM * hmmCon = HMMConventional(2, 2);
    HMM * hmmSparse = HMMConventionalsparse(2, 2);
    HMM * hmmBlas = HMMBLAS(2, 2);
    HMM * hmmCsr = HMMCsr(2, 2);
    
    double transitionProbs[2][2] = {
        {0.5, 0.5},
        {0.3, 0.7}
    };
    
    double emissionProbs[2][2] = {
        {0.3, 0.7},
        {0.8, 0.2}
    };
    
    double initProbs[2] = {0.2, 0.8};
    
    int i;
    int j;
    
    for(i = 0; i < hmmCon->hiddenStates; i++){
        hmmCon->initProbs[i] = initProbs[i];
        hmmBlas->initProbs[i] = initProbs[i];
        hmmCsr->initProbs[i] = initProbs[i];
        hmmSparse->initProbs[i] = initProbs[i];
    }
    for(i = 0; i < hmmCon->hiddenStates; i++){
        for(j = 0; j < hmmCon->hiddenStates; j++){
            hmmCon->transitionProbs[i*hmmCon->hiddenStates+j] = transitionProbs[i][j];
            hmmBlas->transitionProbs[i*hmmCon->hiddenStates+j] = transitionProbs[i][j];
            hmmCsr->transitionProbs[i*hmmCon->hiddenStates+j] = transitionProbs[i][j];
            hmmSparse->transitionProbs[i*hmmCon->hiddenStates+j] = transitionProbs[i][j];
        }
    }
    for(i = 0; i < hmmCon->hiddenStates; i++){
        for(j = 0; j < hmmCon->observations; j++){
            hmmCon->emissionProbs[i*hmmCon->observations+j] = emissionProbs[i][j];
            hmmBlas->emissionProbs[i*hmmCon->observations+j] = emissionProbs[i][j];
            hmmCsr->emissionProbs[i*hmmCon->observations+j] = emissionProbs[i][j];
            hmmSparse->emissionProbs[i*hmmCon->observations+j] = emissionProbs[i][j];
        }
    }
    
    double test[20] = {
        0.838486, 1.015142,
        0.848495, 1.026387,
        0.854859, 1.026767,
        0.898964, 1.018758,
        1.267584, 0.950282,
        1.076550, 0.867237,
        0.944879, 1.143683,
        0.862481, 1.041273,
        0.868282, 1.026152,
        1.000000, 1.000000
    };
    
    const unsigned int observation[10] = {0, 0, 0, 0, 0, 1, 1, 0, 0, 0};
    const unsigned int obsLenght = 10;
    
    double * alpha = calloc(hmmCon->hiddenStates*obsLenght, sizeof(double));
    double * beta = calloc(hmmCon->hiddenStates*obsLenght, sizeof(double));
    double * scaleFactor = calloc(obsLenght, sizeof(double));
    
    F(hmmCon, observation, obsLenght, scaleFactor, alpha);
    B(hmmCon, observation, obsLenght, scaleFactor, beta);
    
    for(i = 0; i < obsLenght; i++){
       for(j = 0; j < hmmCon->hiddenStates; j++){
           assert(fabs(beta[i*hmmCon->hiddenStates+j] - test[i*hmmCon->hiddenStates+j]) < 0.00001);
       }
    }
    
    F(hmmBlas, observation, obsLenght, scaleFactor, alpha);
    B(hmmBlas, observation, obsLenght, scaleFactor, beta);
    
    for(i = 0; i < obsLenght; i++){
       for(j = 0; j < hmmCon->hiddenStates; j++){
           assert(fabs(beta[i*hmmCon->hiddenStates+j] - test[i*hmmCon->hiddenStates+j]) < 0.00001);
       }
    }
    
    F(hmmCsr, observation, obsLenght, scaleFactor, alpha);
    B(hmmCsr, observation, obsLenght, scaleFactor, beta);
    
    for(i = 0; i < obsLenght; i++){
       for(j = 0; j < hmmCon->hiddenStates; j++){
           assert(fabs(beta[i*hmmCon->hiddenStates+j] - test[i*hmmCon->hiddenStates+j]) < 0.00001);
       }
    }
    
    double * scaleFactorsparse = calloc(obsLenght, sizeof(double));
    double * betasparse = calloc(hmmCon->hiddenStates*obsLenght, sizeof(double));

    F(hmmSparse, observation, obsLenght, scaleFactorsparse, alpha);
    B(hmmSparse, observation, obsLenght, scaleFactorsparse, betasparse);
    
    for(i = 0; i < obsLenght; i++){
       for(j = 0; j < hmmCon->hiddenStates; j++){
           assert(fabs(beta[i*hmmCon->hiddenStates+j] - test[i*hmmCon->hiddenStates+j]) < 0.00001);
       }
    }
    
    assert(validateHMM(hmmSparse));
    HMMDeallocate(hmmSparse);
    
    assert(validateHMM(hmmCsr));
    HMMDeallocate(hmmCsr);
    
    assert(validateHMM(hmmCon));
    HMMDeallocate(hmmCon);

    assert(validateHMM(hmmBlas));
    HMMDeallocate(hmmBlas);
    
    free(alpha);
    free(beta);
    free(scaleFactor);
    free(betasparse);
    free(scaleFactorsparse);
    
    return true;
}

bool testBackward7x4(){
    // Conventional and BLAS
    HMM * hmmCon = HMMConventional(7, 4);
    HMM * hmmSparse = HMMConventionalsparse(7, 4);
    HMM * hmmBlas = HMMBLAS(7, 4);
    HMM * hmmCsr = HMMCsr(7, 4);

  double transitionProbs[7][7] = {
          {0.0 , 0.0 , 0.9 , 0.1 , 0.0 , 0.0 , 0.0},
          {1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0},
          {0.0 , 1.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0},
          {0.0 , 0.0 , 0.05 , 0.9 , 0.05 , 0.0 , 0.0},
          {0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0},
          {0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 1.0},
          {0.0 , 0.0 , 0.0 , 0.1 , 0.9 , 0.0 , 0.0}
       };

       double emissionProbs[7][4] = {
          {0.3 , 0.25 , 0.25 , 0.2},
          {0.2 , 0.35 , 0.15 , 0.3},
          {0.4 , 0.15 , 0.2 , 0.25},
          {0.25 , 0.25 , 0.25 , 0.25},
          {0.2 , 0.4 , 0.3 , 0.1},
          {0.3 , 0.2 , 0.3 , 0.2},
          {0.15 , 0.3 , 0.2 , 0.35}
       };

    double initProbs[7] = {0.0 , 0.0 , 0.0 , 1.0 , 0.0 , 0.0 , 0.0};

    int i;
    int j;

    for(i = 0; i < hmmCon->hiddenStates; i++){
        hmmCon->initProbs[i] = initProbs[i];
        hmmBlas->initProbs[i] = initProbs[i];
        hmmCsr->initProbs[i] = initProbs[i];
        hmmSparse->initProbs[i] = initProbs[i];
    }
    for(i = 0; i < hmmCon->hiddenStates; i++){
        for(j = 0; j < hmmCon->hiddenStates; j++){
            hmmCon->transitionProbs[i*hmmCon->hiddenStates+j] = transitionProbs[i][j];
            hmmBlas->transitionProbs[i*hmmCon->hiddenStates+j] = transitionProbs[i][j];
            hmmCsr->transitionProbs[i*hmmCon->hiddenStates+j] = transitionProbs[i][j];
            hmmSparse->transitionProbs[i*hmmCon->hiddenStates+j] = transitionProbs[i][j];
        }
    }
    for(i = 0; i < hmmCon->hiddenStates; i++){
        for(j = 0; j < hmmCon->observations; j++){
            hmmCon->emissionProbs[i*hmmCon->observations+j] = emissionProbs[i][j];
            hmmBlas->emissionProbs[i*hmmCon->observations+j] = emissionProbs[i][j];
            hmmCsr->emissionProbs[i*hmmCon->observations+j] = emissionProbs[i][j];
            hmmSparse->emissionProbs[i*hmmCon->observations+j] = emissionProbs[i][j];
        }
    }

    double test[70] = {1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
        0.631863, 0.987285, 1.382200, 0.997158, 0.789828, 1.184743, 1.520420,
        1.110377, 0.640795, 0.600745, 1.014258, 1.441788, 1.233530, 0.966198,
        0.678220, 0.938276, 0.812214, 1.026372, 1.042341, 1.428777, 0.655376,
        0.536116, 0.671809, 1.301170, 1.021738, 1.132217, 0.779017, 1.588450,
        1.044681, 0.539040, 0.405283, 1.045212, 0.939918, 1.277689, 1.332193,
        0.509689, 0.907717, 0.702552, 1.064127, 1.110176, 2.025684, 0.481032,
        0.733717, 0.404995, 1.081898, 1.008185, 1.609595, 0.668894, 0.502656,
        0.876315, 0.730823, 0.242038, 1.043087, 0.799506, 0.400538, 1.831926,
        0.232683, 0.867638, 1.013022, 1.000000, 0.317258, 2.176545, 1.243165};

    const unsigned int observation[10] = {0,1,2,3,3,2,1,3,2,1};
    const unsigned int obsLenght = 10;

    double * alpha = calloc(hmmCon->hiddenStates*obsLenght, sizeof(double));
    double * beta = calloc(hmmCon->hiddenStates*obsLenght, sizeof(double));
    double * scaleFactor = calloc(obsLenght, sizeof(double));

    F(hmmCon, observation, obsLenght, scaleFactor, alpha);
    B(hmmCon, observation, obsLenght, scaleFactor, beta);

    for(i = 0; i < obsLenght; i++){
       for(j = 0; j < hmmCon->hiddenStates; j++){
           assert(fabs(beta[i*hmmCon->hiddenStates+j] - test[(obsLenght-i-1)*hmmCon->hiddenStates+j]) < 0.00001);
       }
    }

    F(hmmBlas, observation, obsLenght, scaleFactor, alpha);
    B(hmmBlas, observation, obsLenght, scaleFactor, beta);

    for(i = 0; i < obsLenght; i++){
       for(j = 0; j < hmmCon->hiddenStates; j++){
           assert(fabs(beta[i*hmmCon->hiddenStates+j] - test[(obsLenght-i-1)*hmmCon->hiddenStates+j]) < 0.00001);
       }
    }

    F(hmmCsr, observation, obsLenght, scaleFactor, alpha);
    B(hmmCsr, observation, obsLenght, scaleFactor, beta);

    for(i = 0; i < obsLenght; i++){
       for(j = 0; j < hmmCon->hiddenStates; j++){
           assert(fabs(beta[i*hmmCon->hiddenStates+j] - test[(obsLenght-i-1)*hmmCon->hiddenStates+j]) < 0.00001);
       }
    }

    double * scaleFactorsparse = calloc(obsLenght, sizeof(double));
    double * betasparse = calloc(hmmCon->hiddenStates*obsLenght, sizeof(double));

    F(hmmSparse, observation, obsLenght, scaleFactorsparse, alpha);
    B(hmmSparse, observation, obsLenght, scaleFactorsparse, betasparse);

    for(i = 0; i < obsLenght; i++){
       for(j = 0; j < hmmCon->hiddenStates; j++){
           assert(fabs(beta[i*hmmCon->hiddenStates+j] - test[(obsLenght-i-1)*hmmCon->hiddenStates+j]) < 0.00001);
       }
    }

    assert(validateHMM(hmmSparse));
    HMMDeallocate(hmmSparse);

    assert(validateHMM(hmmCsr));
    HMMDeallocate(hmmCsr);

    assert(validateHMM(hmmCon));
    HMMDeallocate(hmmCon);

    assert(validateHMM(hmmBlas));
    HMMDeallocate(hmmBlas);

    free(alpha);
    free(beta);
    free(scaleFactor);
    free(betasparse);
    free(scaleFactorsparse);

    return true;
}

bool testBackwardAlgorithm() {
    
    assert(testBackward2x2());
    assert(testBackward7x4());
    
    return true;
}
