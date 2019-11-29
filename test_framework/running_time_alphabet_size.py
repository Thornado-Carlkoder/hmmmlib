from running_time import *
import binding as hmm_binding
import random, decimal
import time, sys, os


        
print('## Testing varying alphabet size ##', file = sys.stderr)
inputsize = 1500 # the input size is constant. This number will be multiplied with 60 to become relatable with the other tests
start = 2 # the state space
stop = 15 
increment = 2
replicates = 5
file = os.path.dirname(os.path.realpath(__file__)) + "/data/pantro3_X.fasta"
    
        

def standard_test_alphabet(algorithm, hmmType, input_size, start, stop, increment, file= '', algorithm_version = '', **kwargs):
    """ This standard test tests a varying size alphabet with a constant statespace and inputsize."""

    for i in range(start, stop, increment):
        print(f'{algorithm}\t{i}', file = sys.stderr, end = '\t', flush = True)
        test_standard_data = random.choices([j for j in range(i)], k = input_size*60) # Generates a data set with an arbitrary alphabet size (uniform).

        for _ in range(replicates):
            print('r', end = '', file = sys.stderr, flush = True)
            
            o = hmm_binding.binded_HMM(8, i, hmmType = hmmType)
            set_random(o)

            t0 = time.time()
            #test_standard_output = getattr(o, algorithm)(test_standard_data, **kwargs)
            getattr(o, algorithm)(test_standard_data, **kwargs)
            t1 = time.time()
            o.deallocate()

            print(f'alphabetsize, {i}, {t1-t0}, {algorithm}, {o.hmmType}, {algorithm_version}')
        print('', file = sys.stderr, flush = True)

## Conventional ##
standard_test_alphabet("viterbi", "Conventional", inputsize, start, stop, increment, file)
standard_test_alphabet("posteriorDecoding", "Conventional", inputsize, start, stop, increment, file)
standard_test_alphabet("forward", "Conventional", inputsize, start, stop, increment, file)
standard_test_alphabet("backward_time", "Conventional", inputsize, start, stop, increment, file)
standard_test_alphabet("baumWelch", "Conventional", inputsize, start, stop, increment, file, '1', n_iterations = 1)

## Conventional ##
standard_test_alphabet("viterbi", "Conventional sparse", inputsize, start, stop, increment, file)
standard_test_alphabet("posteriorDecoding", "Conventional sparse", inputsize, start, stop, increment, file)
standard_test_alphabet("forward", "Conventional sparse", inputsize, start, stop, increment, file)
standard_test_alphabet("backward_time", "Conventional sparse", inputsize, start, stop, increment, file)
standard_test_alphabet("baumWelch", "Conventional sparse", inputsize, start, stop, increment, file, '1', n_iterations = 1)

## BLAS ##
standard_test_alphabet("viterbi", "BLAS", inputsize, start, stop, increment, file)
standard_test_alphabet("posteriorDecoding", "BLAS", inputsize, start, stop, increment, file)
standard_test_alphabet("forward", "BLAS", inputsize, start, stop, increment, file)
standard_test_alphabet("backward_time", "BLAS", inputsize, start, stop, increment, file)
standard_test_alphabet("baumWelch", "BLAS", inputsize, start, stop, increment, file, '1', n_iterations = 1)

## CSR ##
standard_test_alphabet("viterbi", "CSR", inputsize, start, stop, increment, file)
standard_test_alphabet("posteriorDecoding", "CSR", inputsize, start, stop, increment, file)
standard_test_alphabet("forward", "CSR", inputsize, start, stop, increment, file)
standard_test_alphabet("backward_time", "CSR", inputsize, start, stop, increment, file)
standard_test_alphabet("baumWelch", "CSR", inputsize, start, stop, increment, file, '1', n_iterations = 1)

## RSB ##
standard_test_alphabet("viterbi", "RSB", inputsize, start, stop, increment, file)
standard_test_alphabet("posteriorDecoding", "RSB", inputsize, start, stop, increment, file)
standard_test_alphabet("forward", "RSB", inputsize, start, stop, increment, file)
standard_test_alphabet("backward_time", "RSB", inputsize, start, stop, increment, file)
standard_test_alphabet("baumWelch", "RSB", inputsize, start, stop, increment, file, '1', n_iterations = 1)









        
        
