from running_time import *
import binding as hmm_binding
import random, decimal
import time, sys, os

print('## Testing varying input size ##', file = sys.stderr)
stspace = 7
alphabet = 4
start = 10
stop = 5010 # levels out at 3000
increment = 500
replicates = 5

file = os.path.dirname(os.path.realpath(__file__)) + "/data/pantro3_X.fasta"
#print(file)


def standard_test_inputsize(algorithm, hmmType, stspace, alphabet, start, stop, increment, file, sparseness, algorithm_version = '', linewidth = 60, **kwargs):
    for i in range(start, stop, increment):
        test_standard_data = [i for i in read_fasta(i, file)]
        print(f'{algorithm}\t{i*linewidth}', file = sys.stderr, end = '\t', flush = True)


        o = hmm_binding.binded_HMM(stspace, alphabet, hmmType = hmmType)
        for replicate in range(replicates):
            print('r', end = '', file = sys.stderr, flush = True)    
            
            set_random(o, sparseness)

            t0 = time.time()
            test_standard_output = getattr(o, algorithm)(test_standard_data, **kwargs)
            t1 = time.time()
            
            print(f'inputsize, {i*linewidth}, {t1-t0}, {algorithm}, {o.hmmType}, {algorithm_version}')
        print('', file = sys.stderr, flush = True) # newline
        o.deallocate()





## Conventional #
standard_test_inputsize("viterbi", "Conventional", stspace, alphabet, start, stop, increment, file, 0)
standard_test_inputsize("posteriorDecoding", "Conventional", stspace, alphabet, start, stop, increment, file, 0)
standard_test_inputsize("forward", "Conventional", stspace, alphabet, start, stop, increment, file, 0)
standard_test_inputsize("backward_time", "Conventional", stspace, alphabet, start, stop, increment, file, 0)
standard_test_inputsize("backward_time", "Conventional", stspace, alphabet, start, stop, increment, file, 0)
for i in range(1, 7):
    standard_test_inputsize("baumWelch", "Conventional", stspace, alphabet, start, stop, increment, file, 0, str(i), n_iterations = 1)



## Conventional sparse #
standard_test_inputsize("viterbi", "Conventional sparse", stspace, alphabet, start, stop, increment, file, 0)
standard_test_inputsize("posteriorDecoding", "Conventional sparse", stspace, alphabet, start, stop, increment, file, 0)
standard_test_inputsize("forward", "Conventional sparse", stspace, alphabet, start, stop, increment, file, 0)
standard_test_inputsize("backward_time", "Conventional sparse", stspace, alphabet, start, stop, increment, file, 0)
standard_test_inputsize("backward_time", "Conventional sparse", stspace, alphabet, start, stop, increment, file, 0)
for i in range(1, 7):
    standard_test_inputsize("baumWelch", "Conventional sparse", stspace, alphabet, start, stop, increment, file, 0, str(i), n_iterations = 1)



## BLAS ##
standard_test_inputsize("viterbi", "BLAS", stspace, alphabet, start, stop, increment, file, 0)
standard_test_inputsize("posteriorDecoding", "BLAS", stspace, alphabet, start, stop, increment, file, 0)
standard_test_inputsize("forward", "BLAS", stspace, alphabet, start, stop, increment, file, 0)
standard_test_inputsize("backward_time", "BLAS", stspace, alphabet, start, stop, increment, file, 0)
standard_test_inputsize("backward_time", "BLAS", stspace, alphabet, start, stop, increment, file, 0)
for i in range(1, 7):
    standard_test_inputsize("baumWelch", "BLAS", stspace, alphabet, start, stop, increment, file, 0, str(i), n_iterations = 1)


## CSR ##
standard_test_inputsize("viterbi", "CSR", stspace, alphabet, start, stop, increment, file, 0)
standard_test_inputsize("posteriorDecoding", "CSR", stspace, alphabet, start, stop, increment, file, 0)
standard_test_inputsize("forward", "CSR", stspace, alphabet, start, stop, increment, file, 0)
standard_test_inputsize("backward_time", "CSR", stspace, alphabet, start, stop, increment, file, 0)
standard_test_inputsize("backward_time", "CSR", stspace, alphabet, start, stop, increment, file, 0)
for i in range(1, 7):
    standard_test_inputsize("baumWelch", "CSR", stspace, alphabet, start, stop, increment, file, 0, str(i), n_iterations = 1)


## RSB ##
standard_test_inputsize("viterbi", "RSB", stspace, alphabet, start, stop, increment, file, 0)
standard_test_inputsize("posteriorDecoding", "RSB", stspace, alphabet, start, stop, increment, file, 0)
standard_test_inputsize("forward", "RSB", stspace, alphabet, start, stop, increment, file, 0)
standard_test_inputsize("backward_time", "RSB", stspace, alphabet, start, stop, increment, file, 0)
standard_test_inputsize("backward_time", "RSB", stspace, alphabet, start, stop, increment, file, 0)
for i in range(1, 7):
    standard_test_inputsize("baumWelch", "RSB", stspace, alphabet, start, stop, increment, file, 0, str(i), n_iterations = 1)







    

