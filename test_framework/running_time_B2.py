from running_time import *
import binding as hmm_binding
import random, decimal
import time, sys, os



print('## Testing varying sparseness of transition and emission matrices. ##', file = sys.stderr)
replicates = 3
input_size = 100000
alphabet_size = 4

start = 0
stop = 1.001
increment = 0.1

if True:
    print('test, observations, time, algorithm, variant, states')

for hidden_states in [400, 300, 100, 30, 10]:

    for algorithm in ['forward', 'backward_time']:
        for hmmType in ['Conventional sparse', 'CSR']:
    
            

            # Generate data
            test_standard_data = random.choices([choice for choice in range(alphabet_size)], k = input_size)  # Generates a data set with an arbitrary alphabet_size size (uniform).

            for sparseness in float_range(start, stop, increment): # increasing sparseness
                print(f'{algorithm}\t{round(sparseness,3)}', file = sys.stderr, end = '\t', flush = True)


                for _ in range(replicates):
                    print('r', end = '', file = sys.stderr, flush = True)
                    
                    o = hmm_binding.binded_HMM(hidden_states, alphabet_size, hmmType = hmmType)
                    
                    set_random(o, sparseness)    

                    t0 = time.time()
                    getattr(o, algorithm)(test_standard_data)
                    t1 = time.time()

                    o.deallocate()

                    print(f'sparseness, {round(sparseness, 4)}, {t1-t0}, {algorithm}, {o.hmmType}, {hidden_states}')
                print('', file = sys.stderr, flush = True) # newline
                

