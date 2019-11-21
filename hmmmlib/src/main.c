//
//  main.c
//  hmmmmlib
//
//  Created by Thor Jakobsen on 29/08/2019.
//  Copyright © 2019 Thor Jakobsen. All rights reserved.
//


#include <stdio.h>

#include "TestForward.h"
#include "TestBackward.h"
#include "TestBaumWelch.h"
#include "TestViterbi.h"
#include "TestPosteriorDecoding.h"

int main() {
    
    /*
     Testing the forward algorithm
     */
    if(testForwardAlgorithm()){
        printf("all Forward Algorithm tests succeeded\n");
    }

    /*
     Testing the backward algorithm
    */
    if(testBackwardAlgorithm()){
        printf("all Backward Algorithm tests succeeded\n");
    }
    /*
     Testing the backward algorithm
     */
    if(testBaumWelch()){
        printf("all Baum Welch tests succeeded\n");
    }
    /*
     Testing the viterbi algorithm
     */
    if(testViterbi()){
        printf("all Viterbi tests succeeded\n");
    }

    if(testPosteriorDecoding()){
        printf("all Posterior Decoding tests succeeded\n");
    }
    
    return 0;
}
    
