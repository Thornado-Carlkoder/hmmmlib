from running_time import *
import binding as hmm_binding
import random, decimal
import time, sys, os


        
print('## Testing varying alphabet size ##', file = sys.stderr)

replicates = 5
input_size = 100000 

# Varying alphabet size
start = 4
stop = 45
increment = 10


if True:
    print('test, observations, time, algorithm, variant')



for algorithm in ['posteriorDecoding', 'baumWelch', 'viterbi']:
    for hmmType in ['Conventional', 'Conventional sparse', 'BLAS', 'CSR']: # Don't forget RSB

        #standard_test_alphabet(algorithm, hmmType, input_size, start, stop, increment, file)

        for alphabet_size in range(start, stop, increment): # alphabet size
            print(f'{algorithm}\t{alphabet_size}', file = sys.stderr, end = '\t', flush = True)
            test_standard_data = random.choices([j for j in range(alphabet_size)], k = input_size) # Generates a data set with an arbitrary alphabet size (uniform).

            for _ in range(replicates):
                print('r', end = '', file = sys.stderr, flush = True)
                
                o = hmm_binding.binded_HMM(8, alphabet_size, hmmType = hmmType)
                set_random(o)

                t0 = time.time()
                
                getattr(o, algorithm)(test_standard_data)
                t1 = time.time()
                o.deallocate()

                print(f'alphabetsize, {alphabet_size}, {t1-t0}, {algorithm}, {o.hmmType}')
            print('', file = sys.stderr, flush = True)
