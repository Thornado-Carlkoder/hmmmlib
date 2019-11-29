from running_time import *
import binding as hmm_binding
import random, decimal
import time, sys, os


print('## Testing varying sparseness of transition and emission matrices. ##', file = sys.stderr)
start = 0
stop = 1.001
increment = 0.1
replicates = 5
inputsize = 1500
hidden_states = 16
file = os.path.dirname(os.path.realpath(__file__)) + "/data/pantro3_X.fasta"

def standard_test_sparseness(algorithm, hmmType, inputsize, hidden_states, start, stop, increment, file, algorithm_version = '', **kwargs):
    """ This standard test tests a varying size statespace with a constant alphabet and inputsize."""

    test_standard_data = [i for i in read_fasta(inputsize, file)]
    for i in float_range(start, stop, increment):
        print(f'{algorithm}\t{i}', file = sys.stderr, end = '\t', flush = True)

        for _ in range(replicates):
            print('r', end = '', file = sys.stderr, flush = True)
            
            o = hmm_binding.binded_HMM(hidden_states, 4, hmmType = hmmType)
            
            set_random(o, i)    

            t0 = time.time()
            getattr(o, algorithm)(test_standard_data, **kwargs)
            t1 = time.time()
            o.deallocate()

            print(f'sparseness, {round(i, 4)}, {t1-t0}, {algorithm}, {o.hmmType}, {algorithm_version}')
        print('', file = sys.stderr, flush = True) # newline

## Conventional ##
standard_test_sparseness("viterbi", "Conventional", inputsize, hidden_states, start, stop, increment, file)
standard_test_sparseness("posteriorDecoding", "Conventional", inputsize, hidden_states, start, stop, increment, file)
standard_test_sparseness("forward", "Conventional", inputsize, hidden_states, start, stop, increment, file)
standard_test_sparseness("backward_time", "Conventional", inputsize, hidden_states, start, stop, increment, file)
standard_test_sparseness("baumWelch", "Conventional", inputsize, hidden_states, start, stop, increment, file, str(1), n_iterations = 1)

## Conventional sparse ##
standard_test_sparseness("viterbi", "Conventional sparse", inputsize, hidden_states, start, stop, increment, file)
standard_test_sparseness("posteriorDecoding", "Conventional sparse", inputsize, hidden_states, start, stop, increment, file)
standard_test_sparseness("forward", "Conventional sparse", inputsize, hidden_states, start, stop, increment, file)
standard_test_sparseness("backward_time", "Conventional sparse", inputsize, hidden_states, start, stop, increment, file)
standard_test_sparseness("baumWelch", "Conventional sparse", inputsize, hidden_states, start, stop, increment, file, str(1), n_iterations = 1)


## BLAS #
standard_test_sparseness("viterbi", "BLAS", inputsize, hidden_states, start, stop, increment, file)
standard_test_sparseness("posteriorDecoding", "BLAS", inputsize, hidden_states, start, stop, increment, file)
standard_test_sparseness("forward", "BLAS", inputsize, hidden_states, start, stop, increment, file)
standard_test_sparseness("backward_time", "BLAS", inputsize, hidden_states, start, stop, increment, file)
standard_test_sparseness("baumWelch", "BLAS", inputsize, hidden_states, start, stop, increment, file, str(1), n_iterations = 1)


## CSR #
standard_test_sparseness("viterbi", "CSR", inputsize, hidden_states, start, stop, increment, file)
standard_test_sparseness("posteriorDecoding", "CSR", inputsize, hidden_states, start, stop, increment, file)
standard_test_sparseness("forward", "CSR", inputsize, hidden_states, start, stop, increment, file)
standard_test_sparseness("backward_time", "CSR", inputsize, hidden_states, start, stop, increment, file)
standard_test_sparseness("baumWelch", "CSR", inputsize, hidden_states, start, stop, increment, file, str(1), n_iterations = 1)


## RSB #
standard_test_sparseness("viterbi", "RSB", inputsize, hidden_states, start, stop, increment, file)
standard_test_sparseness("posteriorDecoding", "RSB", inputsize, hidden_states, start, stop, increment, file)
standard_test_sparseness("forward", "RSB", inputsize, hidden_states, start, stop, increment, file)
standard_test_sparseness("backward_time", "RSB", inputsize, hidden_states, start, stop, increment, file)
standard_test_sparseness("baumWelch", "RSB", inputsize, hidden_states, start, stop, increment, file, str(1), n_iterations = 1)








    
        
