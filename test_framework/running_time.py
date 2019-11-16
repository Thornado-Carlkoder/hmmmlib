import binding as hmm_binding
import time, sys, random

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


def random_row(n):
        row = [random.random() for _ in range(n)]
        row_sum = sum(row)
        row = [i/row_sum for i in row]

        return row


def random_matrix(m, n, custom_seed = None):
    """ Generates a random matrix with size m*n. """
    if custom_seed != None:
        random.seed(custom_seed)

    
    return [random_row(n) for _ in range(m)]

def set_random_dense(object):
    """ Sets all the matrices in a given hmm to random values. (Dense) 
        It automatically reads all sizes. """
    object.setInitProbs(random_row(object.n_hiddenstates))
    object.setTransitionProbs(random_matrix(object.n_hiddenstates, object.n_hiddenstates))
    object.setEmissionProbs(random_matrix(object.n_hiddenstates, object.n_observations))
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


def standard_test_inputsize(algorithm, hmmType, stspace, alphabet, start, stop, increment, file= '', algorithm_version = '', linewidth = 60, **kwargs):
    algorithm_name = algorithm#.__name__


    for i in range(start, stop, increment):
        test_standard_data = [i for i in read_fasta(i, file)]
        print(f'{algorithm_name}\t{i*linewidth}', file = sys.stderr, end = '\t', flush = True)


        o = hmm_binding.binded_HMM(stspace, alphabet, hmmType = hmmType)
        for replicate in range(replicates):
            print('r', end = '', file = sys.stderr, flush = True)    
            
            set_random_dense(o)

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
            set_random_dense(o)

            t0 = time.time()
            test_standard_output = getattr(o, algorithm)(test_standard_data, **kwargs)
            t1 = time.time()
            o.deallocate()

            print(f'statespace, {i}, {t1-t0}, {algorithm}, {o.hmmType}, {algorithm_version}')
        print('', file = sys.stderr, flush = True) # newline







if __name__ == "__main__" :
        
    # Make sure that the users specified a test to run.
    arguments = [i for i in sys.argv]

    
    #value_when_true if condition else value_when_false

    #value_when_true if condition else value_when_false

    run_input = True if 'input' in arguments else False
    run_statespace = True if 'statespace' in arguments else False
    run_alphabet = True if 'alphabet' in arguments else False

    if len(arguments) <= 1:
        print("Error: Please give an argument <input|statespace|alphabet>")
        print("example:")
        print("$ python running_time.py input")
    else:
        print('test, observations, time, algorithm, variant, iterations')


    # Test definition for varying input size.
    if run_input:
        print('## Testing varying input size ##', file = sys.stderr)
        stspace = 7
        alphabet = 4
        start = 10
        stop = 2010
        increment = 500
        replicates = 4
        file = '../../test_framework/data/pantro3_X.fasta'


        
        ## Conventional #
        standard_test_inputsize("viterbi", "Conventional", stspace, alphabet, start, stop, increment, file)
        standard_test_inputsize("posteriorDecoding", "Conventional", stspace, alphabet, start, stop, increment, file)
        standard_test_inputsize("forward", "Conventional", stspace, alphabet, start, stop, increment, file)
        standard_test_inputsize("backward", "Conventional", stspace, alphabet, start, stop, increment, file) 
        for i in range(1, 7):
            standard_test_inputsize("baumWelch", "Conventional", stspace, alphabet, start, stop, increment, file, str(i), n_iterations = 1)

        ## BLAS ##
        standard_test_inputsize("viterbi", "BLAS", stspace, alphabet, start, stop, increment, file)
        standard_test_inputsize("posteriorDecoding", "BLAS", stspace, alphabet, start, stop, increment, file)
        standard_test_inputsize("forward", "BLAS", stspace, alphabet, start, stop, increment, file)
        standard_test_inputsize("backward", "BLAS", stspace, alphabet, start, stop, increment, file) 
        for i in range(1, 7):
            standard_test_inputsize("baumWelch", "BLAS", stspace, alphabet, start, stop, increment, file, str(i), n_iterations = 1)

        ## CSR ##
        standard_test_inputsize("viterbi", "CSR", stspace, alphabet, start, stop, increment, file)
        standard_test_inputsize("posteriorDecoding", "CSR", stspace, alphabet, start, stop, increment, file)
        standard_test_inputsize("forward", "CSR", stspace, alphabet, start, stop, increment, file)
        standard_test_inputsize("backward", "CSR", stspace, alphabet, start, stop, increment, file) 
        for i in range(1, 7):
            standard_test_inputsize("baumWelch", "CSR", stspace, alphabet, start, stop, increment, file, str(i), n_iterations = 1)


    # Test definition for varying state space size.
    if run_statespace:
        print('## Testing varying state space ##', file = sys.stderr)
        inputsize = 1000 # the input size is constant
        start = 2 # the state space
        stop = 20
        increment = 4
        replicates = 2
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

        







