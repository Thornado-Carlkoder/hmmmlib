import binding as hmm_binding
import time, sys, random, decimal


def float_range(start, stop, step):
    while start < stop:
        yield float(start)
        start += decimal.Decimal(step)

def read_fasta(n_lines, file):
    mapping = {letter: index for index, letter in enumerate(['A', 'C', 'G', 'T'])}
    
    with open(file, 'r') as fasta_file:
        for line_no, line in enumerate(fasta_file):
            if line[0] == '>':
                continue
            elif n_lines == line_no:
                break
            for number in [mapping[i] for i in str(line).upper().strip()]:
                yield number


def random_row(n, sparseness = 0, seed = None):
        row = [random.uniform(-1, 1) for _ in range(n)] # random numbers uniformly between [-1, 1)

        # sparseness?
        # A sparseness of 1 (sparse) means that only a single value per row is non-zero.
        # A sparseness of 0 (dense) means that all values are non-zero.
        if sparseness > 0:
            chosen = random.sample([p for p in range(n)], int(round((n-1)*sparseness)))    
            for c in chosen:
                row[c] = 0
            


        #if normalize: # Set the rowsums equal to 1
        row_sum = sum(row)
        row = [i/row_sum for i in row]

        return row


def random_matrix(m, n, sparseness = 0, seed = None):
    """ Generates a random matrix with size m*n. """
    if seed != None:
        random.seed(seed)
    
    return [random_row(n, sparseness, seed) for _ in range(m)]






def set_random(object, sparseness = 0):
    """ Sets all the matrices in a given hmm to random values. (Dense) 
        It automatically reads all sizes. """
    object.setInitProbs(random_row(object.n_hiddenstates, sparseness))
    object.setTransitionProbs(random_matrix(object.n_hiddenstates, object.n_hiddenstates, sparseness))
    object.setEmissionProbs(random_matrix(object.n_hiddenstates, object.n_observations, sparseness))
    return

def set_random_sparse(object, percent):
    pass


def set_sparse_1(object):
    object.setInitProbs([0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00])
    object.setTransitionProbs([[0.00, 0.00, 0.90, 0.10, 0.00, 0.00, 0.00],
                        [1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
                        [0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00],
                        [0.00, 0.00, 0.05, 0.90, 0.05, 0.00, 0.00],
                        [0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00],
                        [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00],
                        [0.00, 0.00, 0.00, 0.10, 0.90, 0.00, 0.00]])
    object.setEmissionProbs([[0.30, 0.25, 0.25, 0.20],
                        [0.20, 0.35, 0.15, 0.30],
                        [0.40, 0.15, 0.20, 0.25],
                        [0.25, 0.25, 0.25, 0.25],
                        [0.20, 0.40, 0.30, 0.10],
                        [0.30, 0.20, 0.30, 0.20],
                        [0.15, 0.30, 0.20, 0.35]])


def standard_test_inputsize(algorithm, hmmType, stspace, alphabet, start, stop, increment, file, sparseness, algorithm_version = '', linewidth = 60, **kwargs):
    algorithm_name = algorithm#.__name__


    for i in range(start, stop, increment):
        test_standard_data = [i for i in read_fasta(i, file)]
        print(f'{algorithm_name}\t{i*linewidth}', file = sys.stderr, end = '\t', flush = True)


        o = hmm_binding.binded_HMM(stspace, alphabet, hmmType = hmmType)
        for replicate in range(replicates):
            print('r', end = '', file = sys.stderr, flush = True)    
            
            set_random(o, sparseness)

            t0 = time.time()
            test_standard_output = getattr(o, algorithm)(test_standard_data, **kwargs)
            t1 = time.time()
            
            print(f'inputsize, {i*linewidth}, {t1-t0}, {algorithm_name}, {o.hmmType}, {algorithm_version}')
        print('', file = sys.stderr, flush = True) # newline

        o.deallocate() 


def standard_test_statespace(algorithm, hmmType, inputsize, start, stop, increment, file= '', algorithm_version = '', **kwargs): 
    """ This standard test tests a varying size statespace with a constant alphabet and inputsize."""

    test_standard_data = [i for i in read_fasta(inputsize, file)]
    for i in range(start, stop, increment):
        print(f'{algorithm}\t{i}', file = sys.stderr, end = '\t', flush = True)

        for replicate in range(replicates):
            print('r', end = '', file = sys.stderr, flush = True)    
            
            o = hmm_binding.binded_HMM(i, 4, hmmType = hmmType)
            set_random(o)

            t0 = time.time()
            test_standard_output = getattr(o, algorithm)(test_standard_data, **kwargs)
            t1 = time.time()
            o.deallocate()

            print(f'statespace, {i}, {t1-t0}, {algorithm}, {o.hmmType}, {algorithm_version}')
        print('', file = sys.stderr, flush = True) # newline








if __name__ == "__main__" :


        
    # Make sure that the users specified a test to run.
    given_arguments = [i for i in sys.argv]


    
    if 'input' not in given_arguments \
        and 'statespace' not in given_arguments \
        and 'alphabet' not in given_arguments \
        and 'sparse' not in given_arguments:
        print("Error: Please give an argument <input|statespace|alphabet|sparse>")
        print("example:")
        print("$ python running_time.py input")
    else:
        print('test, observations, time, algorithm, variant, iterations')


    # Test definition for varying input size.
    if 'input' in given_arguments:
        print('## Testing varying input size ##', file = sys.stderr)
        stspace = 7
        alphabet = 4
        start = 10
        stop = 4010 # levels out at 3000
        increment = 500
        replicates = 5
        file = '../../test_framework/data/pantro3_X.fasta'


        
        ## Conventional #
        standard_test_inputsize("viterbi", "Conventional", stspace, alphabet, start, stop, increment, file, 0)
        standard_test_inputsize("posteriorDecoding", "Conventional", stspace, alphabet, start, stop, increment, file, 0)
        standard_test_inputsize("forward", "Conventional", stspace, alphabet, start, stop, increment, file, 0)
        standard_test_inputsize("backward", "Conventional", stspace, alphabet, start, stop, increment, file, 0) 
        for i in range(1, 7):
            standard_test_inputsize("baumWelch", "Conventional", stspace, alphabet, start, stop, increment, file, 0, str(i), n_iterations = 1)

        ## BLAS ##
        standard_test_inputsize("viterbi", "BLAS", stspace, alphabet, start, stop, increment, file, 0)
        standard_test_inputsize("posteriorDecoding", "BLAS", stspace, alphabet, start, stop, increment, file, 0)
        standard_test_inputsize("forward", "BLAS", stspace, alphabet, start, stop, increment, file, 0)
        standard_test_inputsize("backward", "BLAS", stspace, alphabet, start, stop, increment, file, 0) 
        for i in range(1, 7):
            standard_test_inputsize("baumWelch", "BLAS", stspace, alphabet, start, stop, increment, file, 0, str(i), n_iterations = 1)

        ## CSR ##
        standard_test_inputsize("viterbi", "CSR", stspace, alphabet, start, stop, increment, file, 0)
        standard_test_inputsize("posteriorDecoding", "CSR", stspace, alphabet, start, stop, increment, file, 0)
        standard_test_inputsize("forward", "CSR", stspace, alphabet, start, stop, increment, file, 0)
        standard_test_inputsize("backward", "CSR", stspace, alphabet, start, stop, increment, file, 0) 
        for i in range(1, 7):
            standard_test_inputsize("baumWelch", "CSR", stspace, alphabet, start, stop, increment, file, 0, str(i), n_iterations = 1)


    # Test definition for varying state space size.
    if 'statespace' in given_arguments:
        print('## Testing varying state space ##', file = sys.stderr)
        inputsize = 2000 # the input size is constant
        start = 2 # the state space
        stop = 40 # baumwelch cache-jumps at 20
        increment = 4
        replicates = 5
        file = '../../test_framework/data/pantro3_X.fasta'


        
        ## Conventional ##        
        standard_test_statespace("viterbi", "Conventional", inputsize, start, stop, increment, file)        
        standard_test_statespace("posteriorDecoding", "Conventional", inputsize, start, stop, increment, file)
        standard_test_statespace("forward", "Conventional", inputsize, start, stop, increment, file)
        standard_test_statespace("backward", "Conventional", inputsize, start, stop, increment, file) 
        standard_test_statespace("baumWelch", "Conventional", inputsize, start, stop, increment, file, '1', n_iterations = 1)

        ## BLAS ##        
        standard_test_statespace("viterbi", "BLAS", inputsize, start, stop, increment, file)        
        standard_test_statespace("posteriorDecoding", "BLAS", inputsize, start, stop, increment, file)
        standard_test_statespace("forward", "BLAS", inputsize, start, stop, increment, file)
        standard_test_statespace("backward", "BLAS", inputsize, start, stop, increment, file) 
        standard_test_statespace("baumWelch", "BLAS", inputsize, start, stop, increment, file, '1', n_iterations = 1)

        ## CSR ##        
        standard_test_statespace("viterbi", "CSR", inputsize, start, stop, increment, file)        
        standard_test_statespace("posteriorDecoding", "CSR", inputsize, start, stop, increment, file)
        standard_test_statespace("forward", "CSR", inputsize, start, stop, increment, file)
        standard_test_statespace("backward", "CSR", inputsize, start, stop, increment, file) 
        standard_test_statespace("baumWelch", "CSR", inputsize, start, stop, increment, file, '1', n_iterations = 1)




    if 'sparse' in given_arguments:
        def standard_test_sparseness(algorithm, hmmType, inputsize, start, stop, increment, file, algorithm_version = '', **kwargs): 
            """ This standard test tests a varying size statespace with a constant alphabet and inputsize."""

            test_standard_data = [i for i in read_fasta(inputsize, file)]
            for i in float_range(start, stop, increment):
                print(f'{algorithm}\t{i}', file = sys.stderr, end = '\t', flush = True)

                for replicate in range(replicates):
                    print('r', end = '', file = sys.stderr, flush = True)    
                    
                    o = hmm_binding.binded_HMM(7, 4, hmmType = hmmType)
                    set_random(o, i)

                    t0 = time.time()
                    test_standard_output = getattr(o, algorithm)(test_standard_data, **kwargs)
                    t1 = time.time()
                    o.deallocate()

                    print(f'sparseness, {round(i, 4)}, {t1-t0}, {algorithm}, {o.hmmType}, {algorithm_version}')
                print('', file = sys.stderr, flush = True) # newline






        print('## Testing varying sparseness of transition and emission matrices. ##', file = sys.stderr)
        start = 0
        stop = 1.1
        increment = 0.2
        replicates = 2
        inputsize = 1000
        file = '../../test_framework/data/pantro3_X.fasta'

        
        ## Conventional #
        standard_test_sparseness("viterbi", "Conventional", inputsize, start, stop, increment, file)
        standard_test_sparseness("posteriorDecoding", "Conventional", inputsize, start, stop, increment, file)
        standard_test_sparseness("forward", "Conventional", inputsize, start, stop, increment, file)
        standard_test_sparseness("backward", "Conventional", inputsize, start, stop, increment, file) 
        for i in range(1, 7):
            standard_test_sparseness("baumWelch", "Conventional", inputsize, start, stop, increment, file, str(i), n_iterations = 1)
        