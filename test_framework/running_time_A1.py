from running_time import *
import binding as hmm_binding
import random, decimal
import time, sys, os

print('Testing varying input size', file = sys.stderr)

alphabet_size = 4
replicates = 5
start = 10
stop = 1000010  # levels out at ~ 200K
increment = 100000

if True: # print csv-header
    print('test, observations, time, algorithm, variant, statespace')

for stspace in [300, 100, 30, 10]:
    for algorithm in ['forward', 'backward_time']:
        for hmmType in ['Conventional', 'BLAS']:
            for input_size in range(start, stop, increment):
            
                # progress to stdout
                print(f'{algorithm}\t{input_size}', file = sys.stderr, end = '\t', flush = True)
                
                # Generate data
                test_standard_data = random.choices([choice for choice in range(alphabet_size)], k = input_size)  # Generates a data set with an arbitrary alphabet_size size (uniform).

                o = hmm_binding.binded_HMM(stspace, alphabet_size, hmmType = hmmType)
                for _ in range(replicates):
                    print('r', end = '', file = sys.stderr, flush = True)    
                    
                    set_random(o, sparseness = 0)

                    t0 = time.time()
                    test_standard_output = getattr(o, algorithm)(test_standard_data)
                    t1 = time.time()
                    
                    # Print to csv
                    print(f'inputsize, {input_size}, {t1-t0}, {algorithm}, {o.hmmType}, {""}')
                
                print('', file = sys.stderr, flush = True) # newline
                o.deallocate()