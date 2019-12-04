from running_time import *
import binding as hmm_binding
import random, decimal
import time, sys, os

## Varying state space size ##
print('## Testing varying state space ##', file = sys.stderr)
inputsize = 2000 # the input size is constant
start = 2 # the state space
stop = 40 # baumwelch cache-jumps at 20
increment = 4
replicates = 5
file = os.path.dirname(os.path.realpath(__file__)) + "/data/pantro3_X.fasta"


def standard_test_statespace(algorithm, hmmType, inputsize, start, stop, increment, file= '', algorithm_version = '', **kwargs):
    """ This standard test tests a varying size statespace with a constant alphabet and inputsize."""

    test_standard_data = [i for i in read_fasta(inputsize, file)]
    for i in range(start, stop, increment):
        print(f'{algorithm}\t{i}', file = sys.stderr, end = '\t', flush = True)

        for _ in range(replicates):
            print('r', end = '', file = sys.stderr, flush = True)    
            
            o = hmm_binding.binded_HMM(i, 4, hmmType = hmmType)
            set_random(o)

            t0 = time.time()
            getattr(o, algorithm)(test_standard_data, **kwargs)
            t1 = time.time()
            o.deallocate()

            print(f'statespace, {i}, {t1-t0}, {algorithm}, {o.hmmType}, {algorithm_version}')
        print('', file = sys.stderr, flush = True) # newline





## Conventional ##        
standard_test_statespace("viterbi", "Conventional", inputsize, start, stop, increment, file)
standard_test_statespace("posteriorDecoding", "Conventional", inputsize, start, stop, increment, file)
standard_test_statespace("forward", "Conventional", inputsize, start, stop, increment, file)
standard_test_statespace("backward_time", "Conventional", inputsize, start, stop, increment, file)
standard_test_statespace("baumWelch", "Conventional", inputsize, start, stop, increment, file, '1', n_iterations = 1)

## Conventional sparse ##        
standard_test_statespace("viterbi", "Conventional sparse", inputsize, start, stop, increment, file)
standard_test_statespace("posteriorDecoding", "Conventional sparse", inputsize, start, stop, increment, file)
standard_test_statespace("forward", "Conventional sparse", inputsize, start, stop, increment, file)
standard_test_statespace("backward_time", "Conventional sparse", inputsize, start, stop, increment, file)
standard_test_statespace("baumWelch", "Conventional sparse", inputsize, start, stop, increment, file, '1', n_iterations = 1)

## BLAS ##        
standard_test_statespace("viterbi", "BLAS", inputsize, start, stop, increment, file)
standard_test_statespace("posteriorDecoding", "BLAS", inputsize, start, stop, increment, file)
standard_test_statespace("forward", "BLAS", inputsize, start, stop, increment, file)
standard_test_statespace("backward_time", "BLAS", inputsize, start, stop, increment, file)
standard_test_statespace("baumWelch", "BLAS", inputsize, start, stop, increment, file, '1', n_iterations = 1)

## CSR ##        
standard_test_statespace("viterbi", "CSR", inputsize, start, stop, increment, file)
standard_test_statespace("posteriorDecoding", "CSR", inputsize, start, stop, increment, file)
standard_test_statespace("forward", "CSR", inputsize, start, stop, increment, file)
standard_test_statespace("backward_time", "CSR", inputsize, start, stop, increment, file)
standard_test_statespace("baumWelch", "CSR", inputsize, start, stop, increment, file, '1', n_iterations = 1)

## RSB ##        
standard_test_statespace("viterbi", "RSB", inputsize, start, stop, increment, file)
standard_test_statespace("posteriorDecoding", "RSB", inputsize, start, stop, increment, file)
standard_test_statespace("forward", "RSB", inputsize, start, stop, increment, file)
standard_test_statespace("backward_time", "RSB", inputsize, start, stop, increment, file)
standard_test_statespace("baumWelch", "RSB", inputsize, start, stop, increment, file, '1', n_iterations = 1)






