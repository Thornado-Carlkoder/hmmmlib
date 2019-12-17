from running_time import *
import binding as hmm_binding
import random, decimal
import time, sys


"""
Figur C 4: Baum-Welch iterations VS scaled time

Why: Sanity check; does increasing number of iterations increase by a constant?
"""



input_size = 100000 # the input size is constant
alphabet_size = 4
replicates = 3
hidden_states = 10

# Number of iterations
start = 2
stop = 20
increment = 1



if True:
    print('test, observations, time, algorithm, variant')
    

for hmmType in ['Conventional sparse', 'Conventional', 'BLAS', 'CSR', 'RSB']: # Don't forget RSB
    
    test_standard_data = random.choices([j for j in range(alphabet_size)], k = input_size) # Generates a data set with an arbitrary alphabet size (uniform).
    
    for n_iterations in range(start, stop, increment):
        print(f'{"baumWelch"}\t{n_iterations}', file = sys.stderr, end = '\t', flush = True)

        for _ in range(replicates):
            print('r', end = '', file = sys.stderr, flush = True)    

            o = hmm_binding.binded_HMM(n_iterations, 4, hmmType = hmmType)
            set_random(o)

            t0 = time.time()
            getattr(o, 'baumWelch')(test_standard_data, n_iterations = n_iterations)
            t1 = time.time()
            o.deallocate()

            print(f'baumwelchiterations, {n_iterations}, {t1-t0}, {"baumWelch"}, {o.hmmType}')
        print('', file = sys.stderr, flush = True) # newline




