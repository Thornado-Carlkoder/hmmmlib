#include "TestBackward.h"
#include "hmm.h"

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

extern bool testBackwardAlgorithm() {
    HMM * hmmCon = HMMConventional(2, 2);
    HMM * hmmBlas = HMMConventional(2, 2);
    
    double transitionProbs[2][2] = {
        {0.5, 0.5},
        {0.3, 0.7}
    };
    
    double emissionProbs[2][2] = {
        {0.3, 0.7},
        {0.8, 0.2}
    };
    
    double initProbs[2] = {0.2, 0.8};
    
    for(unsigned int i = 0; i < hmmCon->hiddenStates; i++){
        hmmCon->initProbs[i] = initProbs[i];
        hmmBlas->initProbs[i] = initProbs[i];
    }
    int i;
    int j;
    for(i = 0; i < hmmCon->hiddenStates; i++){
        for(j = 0; j < hmmCon->hiddenStates; j++){
            hmmCon->transitionProbs[i*hmmCon->hiddenStates+j] = transitionProbs[i][j];
            hmmBlas->transitionProbs[i*hmmCon->hiddenStates+j] = transitionProbs[i][j];
        }
    }
    for(i = 0; i < hmmCon->hiddenStates; i++){
        for(j = 0; j < hmmCon->observations; j++){
            hmmCon->emissionProbs[i*hmmCon->observations+j] = emissionProbs[i][j];
            hmmBlas->emissionProbs[i*hmmCon->observations+j] = emissionProbs[i][j];
        }
    }
    
    const unsigned int observation[10] = {0, 0, 0, 0, 0, 1, 1, 0, 0, 0};
    const unsigned int obsLenght = 10;
    double * alphaCon = calloc(hmmCon->hiddenStates*obsLenght, sizeof(double));
    double * betaCon = calloc(hmmCon->hiddenStates*obsLenght, sizeof(double));
    
    double * alphaBlas = calloc(hmmCon->hiddenStates*obsLenght, sizeof(double));
    double * betaBlas = calloc(hmmCon->hiddenStates*obsLenght, sizeof(double));
    
    double * scaleFactorAlpha = calloc(obsLenght, sizeof(double));
    F(hmmCon, observation, obsLenght, scaleFactorAlpha, alphaCon);
    B(hmmCon, observation, obsLenght, scaleFactorAlpha, betaCon);
    
    double * scaleFactorBlas = calloc(obsLenght, sizeof(double));
    F(hmmBlas, observation, obsLenght, scaleFactorBlas, alphaBlas);
    B(hmmBlas, observation, obsLenght, scaleFactorBlas, betaBlas);
    
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
    
    
    for(i = 0; i < obsLenght; i++){
       for(j = 0; j < hmmCon->hiddenStates; j++){
           assert(fabs(betaBlas[i*hmmCon->hiddenStates+j] - test[i*hmmCon->hiddenStates+j]) < 0.00001);
           assert(fabs(betaCon[i*hmmCon->hiddenStates+j] - test[i*hmmCon->hiddenStates+j]) < 0.00001);
       }
    }
    
    assert(validateHMM(hmmCon) == true);
    HMMDeallocate(hmmCon);

    assert(validateHMM(hmmBlas) == true);
    HMMDeallocate(hmmBlas);
    
    HMM * hmm2 = HMMBLAS(7, 4);
           
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

    const unsigned int observation2[10] = {0,1,2,3,3,2,1,3,2,1};
    const unsigned int obsLenght2 = 10;

    double * scaleFactor2 = calloc(obsLenght2, sizeof(double));
    double * alpha2 = calloc(obsLenght2*hmm2->hiddenStates, sizeof(double));
    F(hmm2, observation2, obsLenght2, scaleFactor2, alpha2);
    double * beta2 = calloc(obsLenght2*hmm2->hiddenStates, sizeof(double));
    B(hmm2, observation2, obsLenght2, scaleFactor2, beta2);
    
    double test2alpha[70] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0297029702970297, 0.8910891089108911, 0.07920792079207922, 0.0, 0.0,
        0.0, 0.017751479289940822, 0.03550295857988165, 0.7988165680473371, 0.05325443786982248, 0.09467455621301775, 0.0,
        0.014105201292976783, 0.04231560387893035, 0.03967087863649721, 0.7140758154569498, 0.015868351454598886, 0.04231560387893036, 0.13164854540111665,
        0.0367677567200097, 0.05170465788751365, 0.05256640218563887, 0.7138434434025496, 0.06698625677426767, 0.013787908770003642, 0.06434357426001698,
        0.05198659387687777, 0.03171182226489543, 0.055326572081016934, 0.6561285776082829, 0.112934136194607, 0.08082182461391929, 0.011090473360400593,
        0.03141205283200857, 0.07672500391140769, 0.04730517836050952, 0.5911817003208785, 0.06781341533059485, 0.08949326274723163, 0.09606938649736914,
        0.06483316415628948, 0.05995980914443716, 0.06108337113044319, 0.575462087135659, 0.049019502784910594, 0.05730280956653559, 0.13233925608172487,
        0.06080739819354018, 0.03716810567944167, 0.07068361444879583, 0.5452330770139235, 0.17996260493700922, 0.059654928211884464, 0.046490271515404996,
        0.03669552875343646, 0.09769886257539809, 0.048567519065541305, 0.49506394220074496, 0.10915885400232408, 0.1421395645172638, 0.07067572888529124};
    
    
    double test2beta[70] = {1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000, 1.000000,
        0.631863, 0.987285, 1.382200, 0.997158, 0.789828, 1.184743, 1.520420,
        1.110377, 0.640795, 0.600745, 1.014258, 1.441788, 1.233530, 0.966198,
        0.678220, 0.938276, 0.812214, 1.026372, 1.042341, 1.428777, 0.655376,
        0.536116, 0.671809, 1.301170, 1.021738, 1.132217, 0.779017, 1.588450,
        1.044681, 0.539040, 0.405283, 1.045212, 0.939918, 1.277689, 1.332193,
        0.509689, 0.907717, 0.702552, 1.064127, 1.110176, 2.025684, 0.481032,
        0.733717, 0.404995, 1.081898, 1.008185, 1.609595, 0.668894, 0.502656,
        0.876315, 0.730823, 0.242038, 1.043087, 0.799506, 0.400538, 1.831926,
        0.232683, 0.867638, 1.013022, 1.000000, 0.317258, 2.176545, 1.243165};
    

    for(i = 0; i < obsLenght2; i++){
        for(j = 0; j < hmm2->hiddenStates; j++){
            assert(fabs(beta2[i*hmm2->hiddenStates+j]-test2beta[(obsLenght2-i-1)*hmm2->hiddenStates+j]) < 0.00001);
        }
    }
    
    for (j = 0; j < obsLenght2; j++){
        for(i = 0; i < hmm2->hiddenStates; i++) {
            assert(fabs(alpha2[j*hmm2->hiddenStates+i]-test2alpha[j*hmm2->hiddenStates+i]) < 0.00001);
        }
    }
    free(alpha2);
    free(beta2);
    free(scaleFactor2);
    assert(validateHMM(hmm2) == true);
    HMMDeallocate(hmm2);
    
    HMM * hmm3 = HMMCsr(7, 4);
    
    for(i = 0; i < hmm3->hiddenStates; i++){
       hmm3->initProbs[i] = initProbs2[i];
    }
    for(i = 0; i < hmm3->hiddenStates; i++){
      for(j = 0; j < hmm3->hiddenStates; j++){
          hmm3->transitionProbs[i*hmm3->hiddenStates+j] = transitionProbs2[i][j];
      }
    }
    for(i = 0; i < hmm3->hiddenStates; i++){
      for(j = 0; j < hmm3->observations; j++){
          hmm3->emissionProbs[i*hmm3->observations+j] = emissionProbs2[i][j];
      }
    }

    double * scaleFactor3 = calloc(obsLenght2, sizeof(double));
    double * alpha3 = calloc(obsLenght2*hmm3->hiddenStates, sizeof(double));
    F(hmm3, observation2, obsLenght2, scaleFactor3, alpha3);
    double * beta3 = calloc(obsLenght2*hmm3->hiddenStates, sizeof(double));
    B(hmm3, observation2, obsLenght2, scaleFactor3, beta3);
    
    for(i = 0; i < obsLenght2; i++){
        for(j = 0; j < hmm2->hiddenStates; j++){
            assert(fabs(beta3[i*hmm2->hiddenStates+j]-test2beta[(obsLenght2-i-1)*hmm2->hiddenStates+j]) < 0.00001);
        }
    }

    assert(validateHMM(hmm3) == true);
    HMMDeallocate(hmm3);
    free(alpha3);
    free(scaleFactor3);
    free(beta3);

    return true;
}
