from running_time import *
import binding as hmm_binding
import random, decimal
import time, sys, os


"""
Figurer  C 3: All algorithms: Varying statespace

Why: Sanity check; does increasing state space increase running time by the square?
"""


print('## Testing varying state space ##', file = sys.stderr)
input_size = 100000 # the input size is constant
alphabet_size = 4

start = 2 # the state space
stop = 100 # baumwelch cache-jumps at 20
increment = 4
replicates = 3



if True:
    print('test, observations, time, algorithm, variant')
    

for algorithm in ['forward', 'backward_time', 'posteriorDecoding', 'baumWelch', 'viterbi']:
    for hmmType in ['Conventional sparse', 'Conventional', 'BLAS', 'CSR', 'RSB']: # Don't forget RSB
        
        test_standard_data = random.choices([j for j in range(alphabet_size)], k = input_size) # Generates a data set with an arbitrary alphabet size (uniform).
        
        for hidden_states in range(start, stop, increment):
            print(f'{algorithm}\t{hidden_states}', file = sys.stderr, end = '\t', flush = True)

            for _ in range(replicates):
                print('r', end = '', file = sys.stderr, flush = True)    
                
                o = hmm_binding.binded_HMM(hidden_states, 4, hmmType = hmmType)
                set_random(o)

                t0 = time.time()
                getattr(o, algorithm)(test_standard_data)
                t1 = time.time()
                o.deallocate()

                print(f'statespace, {hidden_states}, {t1-t0}, {algorithm}, {o.hmmType}')
            print('', file = sys.stderr, flush = True) # newline




