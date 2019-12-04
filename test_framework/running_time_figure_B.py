from running_time import *
import binding as hmm_binding
import random, decimal
import time, sys, os


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

            print(f'sparseness_{inputsize}, {round(i, 4)}, {t1-t0}, {algorithm}, {o.hmmType}, {algorithm_version}')
        print('', file = sys.stderr, flush = True) # newline


print('## Testing varying sparseness of transition and emission matrices. ##', file = sys.stderr)
start = 0
stop = 1.001
increment = 0.1
replicates = 1
inputsize = 1500
file = os.path.dirname(os.path.realpath(__file__)) + "/data/pantro3_X.fasta"


for hidden_states in [400]:#300, 100, 30, 10]:

    print('# hs:', hidden_states, file = sys.stderr)
    
    # Conventional #
    #standard_test_sparseness("forward", "Conventional", inputsize, hidden_states, start, stop, increment, file, hidden_states)
    #standard_test_sparseness("backward_time", "Conventional", inputsize, hidden_states, start, stop, increment, file, hidden_states)

    ## Conventionalsparse #
    #standard_test_sparseness("forward", "Conventional sparse", inputsize, hidden_states, start, stop, increment, file, hidden_states)
    #standard_test_sparseness("backward_time", "Conventional sparse", inputsize, hidden_states, start, stop, increment, file, hidden_states)


    # BLAS #
    standard_test_sparseness("forward", "BLAS", inputsize, hidden_states, start, stop, increment, file, hidden_states)
    standard_test_sparseness("backward_time", "BLAS", inputsize, hidden_states, start, stop, increment, file, hidden_states)

    # CSR ##
    standard_test_sparseness("forward", "CSR", inputsize, hidden_states, start, stop, increment, file, hidden_states)
    standard_test_sparseness("backward_time", "CSR", inputsize, hidden_states, start, stop, increment, file, hidden_states)

    # RSB ##
    standard_test_sparseness("forward", "RSB", inputsize, hidden_states, start, stop, increment, file, hidden_states)
    standard_test_sparseness("backward_time", "RSB", inputsize, hidden_states, start, stop, increment, file, hidden_states)

    
    
    




