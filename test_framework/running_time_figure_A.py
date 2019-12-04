from running_time import *
import binding as hmm_binding
import random, decimal
import time, sys, os

print('Testing varying input size for figure A', file = sys.stderr)
stspace = 8
alphabet = 4
start = 10
stop = 3010 # levels out at 3000
increment = 500
replicates = 5

file = os.path.dirname(os.path.realpath(__file__)) + "/data/pantro3_X.fasta"
#print(file)


def standard_test_inputsize(algorithm, hmmType, stspace, alphabet, start, stop, increment, file, sparseness, algorithm_version = '', linewidth = 60, **kwargs):
    for i in range(start, stop, increment):
        test_standard_data = [i for i in read_fasta(i, file)]
        print(f'{algorithm}\t{i*linewidth}', file = sys.stderr, end = '\t', flush = True)


        o = hmm_binding.binded_HMM(stspace, alphabet, hmmType = hmmType)
        for _ in range(replicates):
            print('r', end = '', file = sys.stderr, flush = True)    
            
            set_random(o, sparseness)

            t0 = time.time()
            test_standard_output = getattr(o, algorithm)(test_standard_data, **kwargs)
            t1 = time.time()
            
            print(f'inputsize, {i*linewidth}, {t1-t0}, {algorithm}, {o.hmmType}, {algorithm_version}, {stspace}')
        print('', file = sys.stderr, flush = True) # newline
        o.deallocate()




for stspace in [16, 8]:
    ## Conventional #
    #standard_test_inputsize("viterbi", "Conventional", stspace, alphabet, start, stop, increment, file, 0)
    standard_test_inputsize("posteriorDecoding", "Conventional", stspace, alphabet, start, stop, increment, file, 0)
    standard_test_inputsize("forward", "Conventional", stspace, alphabet, start, stop, increment, file, 0)
    standard_test_inputsize("backward_time", "Conventional", stspace, alphabet, start, stop, increment, file, 0)
    standard_test_inputsize("backward_time", "Conventional", stspace, alphabet, start, stop, increment, file, 0)
    for i in range(1, 7):
        standard_test_inputsize("baumWelch", "Conventional", stspace, alphabet, start, stop, increment, file, 0, str(i), n_iterations = 1)



    ## Conventional sparse #
    #standard_test_inputsize("viterbi", "Conventional sparse", stspace, alphabet, start, stop, increment, file, 0)
    standard_test_inputsize("posteriorDecoding", "Conventional sparse", stspace, alphabet, start, stop, increment, file, 0)
    standard_test_inputsize("forward", "Conventional sparse", stspace, alphabet, start, stop, increment, file, 0)
    standard_test_inputsize("backward_time", "Conventional sparse", stspace, alphabet, start, stop, increment, file, 0)
    standard_test_inputsize("backward_time", "Conventional sparse", stspace, alphabet, start, stop, increment, file, 0)
    for i in range(1, 7):
        standard_test_inputsize("baumWelch", "Conventional sparse", stspace, alphabet, start, stop, increment, file, 0, str(i), n_iterations = 1)



    ## BLAS ##
    #standard_test_inputsize("viterbi", "BLAS", stspace, alphabet, start, stop, increment, file, 0)
    standard_test_inputsize("posteriorDecoding", "BLAS", stspace, alphabet, start, stop, increment, file, 0)
    standard_test_inputsize("forward", "BLAS", stspace, alphabet, start, stop, increment, file, 0)
    standard_test_inputsize("backward_time", "BLAS", stspace, alphabet, start, stop, increment, file, 0)
    standard_test_inputsize("backward_time", "BLAS", stspace, alphabet, start, stop, increment, file, 0)
    for i in range(1, 7):
        standard_test_inputsize("baumWelch", "BLAS", stspace, alphabet, start, stop, increment, file, 0, str(i), n_iterations = 1)




    """ 
    for algorithm in ['posteriorDecoding', 'forward', 'backward_time', 'backward_time', 'baumWelch']:
        for implementation in ['Conventional', 'Conventional sparse', 'BLAS']:
            for stspace in [10, 30]:
                
                standard_test_inputsize(algorithm, implementation, stspace, alphabet, start, stop, increment, file, 0)
    """




    

