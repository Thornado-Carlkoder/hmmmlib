//
//  main.c
//  hmmmmlib
//
//  Created by Thor Jakobsen on 29/08/2019.
//  Copyright © 2019 Thor Jakobsen. All rights reserved.
//

#include <stdio.h>

#include "hmm.h"
#include "hmm.c"

// Conventional implementations
#include "forward.h"
#include "forward.c"

#include "backward.h"
#include "backward.c"

#include "viterbi.h"
#include "viterbi.c"

#include "baumWelch.h"
#include "baumWelch.c"

#include "posteriorDecoding.h"
#include "posteriorDecoding.c"

// Conventional sparse implementations
#include "forward_con_sparse.c"
#include "forward_con_sparse.h"

#include "backward_con_sparse.c"
#include "backward_con_sparse.h"

// BLAS implementations
#include "forward_blas.h"
#include "forward_blas.c"

#include "backward_blas.h"
#include "backward_blas.c"

// CSR implementations
#include "forward_csr.c"
#include "forward_csr.h"

#include "backward_csr.c"
#include "backward_csr.h"

// RSB (sblas) implementations
#include "forward_sblas.c"
#include "forward_sblas.h"

#include "backward_sblas.c"
#include "backward_sblas.h"



//#include <cblas.h>


