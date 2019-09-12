#include "TestViterbi.h"
#include "viterbi.h"
#include <stdlib.h>
#include <assert.h>

bool testViterbi(){
     HMM * hmm = HMMCreate(2, 3);
    
    
    // set hardcoded probs
    int n_states = 2; // healthy, fever
    int n_obs = 3; // normal, cold, dizzy
    
    // initProbs
    double start_p_[2] = {0.6, 0.4}; 
    
    // transitionProbs
    double trans_p_[2][2] = {{0.7, 0.3},
                             {0.4, 0.6}};
    double** trans_p = calloc(2, sizeof(double*));
    for (int i = 0; i < 2; i++) trans_p[i] = trans_p_[i]; // Hvornår allokerer vi kolonnerne?
    
    // emissionProbs
    double emit_p_[2][3] = {{0.5, 0.4, 0.1},
                            {0.1, 0.3, 0.6}}; // [state][obs]
    double** emit_p = calloc(2, sizeof(double*));
    for (int i = 0; i < 2; i++) emit_p[i] = emit_p_[i];
    
    
    // feed data to the viterbi alg.
    int n_data = 7;
    int data_[7] = {0, 1, 2, 2, 1, 1, 0}; // normal, cold, dizzy, dizzy, cold, cold, normal
    int* data = data_;
    
    int output_[7];
    
   
    
    double** rv = viterbi(hmm, n_obs, n_states, (double*)&start_p_, trans_p, emit_p, n_data, data, (int*)&output_);
    if (rv[0][6] == 0.000035) {
        printf("men det passer jo");
    }
    printf("WGC %f\n", rv[0][6]);
    assert((double)rv[0][6] ==  (double)0.000035);

    


    printf("Viterbi decoding gives the following state sequence:\n");
    for (int i = 0; i < n_data; i++)
    {
        printf("%i, ", output_[i]);
    }
    printf("\n");
    
    for (int i = 0; i < n_data; ++i)
    {
        printf("%f %f \n", rv[0][i], rv [1][i]);
    }
    
    
    
    return true;
}