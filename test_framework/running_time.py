import binding as hmm_binding
import random, decimal
import time, sys

"""
    Description:

    This file imports the defined python binding (hmm_binding) and calls the algorithms with
    varying parameters and records the time it takes for each function/algorithm call.

    This file should be called in the hmmlib/build/ directory with the following command:
        $ python ../../test_framework/running_time.py <input|statespace|sparse|alphabet> | cat

        example:
        $ python ../../test_framework/running_time.py input | cat
        ... which will run the varying input size tests

    Note that this program will output to stderr and stdout at the same time. As a consequence outputs must be piped correctly:
    You can either use the `| cat` or you can pipe it directly to a .csv file for further processing in plotting software.

        example:
        $ python ../../test_framework/running_time.py input > inputsize_runningtime_data_for_future_plotting.csv


    
    Important note: as backward is dependent on the scalefactor from forward, it is very important to always run forward before running backward.


"""

def float_range(start, stop, step):
    """ This function is the float counterpart to range().
        range() will only operate on integers, whereas this 
        function will work on floats. """
    while start < stop:
        yield float(start)
        start += decimal.Decimal(step)


def read_fasta(n_lines, file):
    """ Read n_lines number of lines from a faste file file. 
        Yields one letter for each iteration. """
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
    """ Generates a random row of size n with sparseness between 0 and 1.
        A sparseness of 1 (sparse) means that only a single value per row is non-zero.
        A sparseness of 0 (dense) means that all values are non-zero. """

    row = [random.uniform(0, 1) for _ in range(n)] # random numbers uniformly between [-1, 1)
    
    if sparseness > 0:
        chosen = random.sample([p for p in range(n)], int(round((n-1)*sparseness)))
        for c in chosen:
            row[c] = 0

    row_sum = sum(row)
    row = [i/row_sum for i in row]

    return row


def random_matrix(m, n, sparseness = 0, seed = None):
    """ Generates a random matrix with size m*n. """
    if seed is not None:
        random.seed(seed)
    
    return [random_row(n, sparseness, seed) for _ in range(m)]


def set_random(hmm_obj, sparseness = 0):
    """ Sets all the matrices in a given hmm to random values. (Dense)
        It automatically reads all sizes. """
    hmm_obj.setInitProbs(random_row(hmm_obj.n_hiddenstates, sparseness))
    hmm_obj.setTransitionProbs(random_matrix(hmm_obj.n_hiddenstates, hmm_obj.n_hiddenstates, sparseness))
    hmm_obj.setEmissionProbs(random_matrix(hmm_obj.n_hiddenstates, hmm_obj.n_observations, sparseness))
    return


def set_sparse_1(hmm_obj):
    """ This is an arbitrary sparse set up. """
    hmm_obj.setInitProbs([0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00])
    hmm_obj.setTransitionProbs([[0.00, 0.00, 0.90, 0.10, 0.00, 0.00, 0.00],
                        [1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
                        [0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00],
                        [0.00, 0.00, 0.05, 0.90, 0.05, 0.00, 0.00],
                        [0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00],
                        [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00],
                        [0.00, 0.00, 0.00, 0.10, 0.90, 0.00, 0.00]])
    hmm_obj.setEmissionProbs([[0.30, 0.25, 0.25, 0.20],
                        [0.20, 0.35, 0.15, 0.30],
                        [0.40, 0.15, 0.20, 0.25],
                        [0.25, 0.25, 0.25, 0.25],
                        [0.20, 0.40, 0.30, 0.10],
                        [0.30, 0.20, 0.30, 0.20],
                        [0.15, 0.30, 0.20, 0.35]])




if __name__ == "__main__" :
        
    # Make sure that the users specified a test to run.
    given_arguments = [i for i in sys.argv]
    if 'input' not in given_arguments \
        and 'statespace' not in given_arguments \
        and 'sparse' not in given_arguments \
        and 'alphabet' not in given_arguments:
        print("Error: Please give an argument <input|statespace|alphabet|sparse>")
        print("example:")
        print("$ python running_time.py input")
    else:
        pass
        #print('test, observations, time, algorithm, variant, iterations')


    ## Varying input size ##
    if 'input' in given_arguments:

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

        ## RSB ##
        standard_test_inputsize("viterbi", "RSB", stspace, alphabet, start, stop, increment, file, 0)
        standard_test_inputsize("posteriorDecoding", "RSB", stspace, alphabet, start, stop, increment, file, 0)
        standard_test_inputsize("forward", "RSB", stspace, alphabet, start, stop, increment, file, 0)
        standard_test_inputsize("backward", "RSB", stspace, alphabet, start, stop, increment, file, 0) 
        for i in range(1, 7):
            standard_test_inputsize("baumWelch", "RSB", stspace, alphabet, start, stop, increment, file, 0, str(i), n_iterations = 1)

        

    ## Varying state space size ##
    if 'statespace' in given_arguments:

        def standard_test_statespace(algorithm, hmmType, inputsize, start, stop, increment, file= '', algorithm_version = '', **kwargs):
            """ This standard test tests a varying size statespace with a constant alphabet and inputsize."""

            test_standard_data = [i for i in read_fasta(inputsize, file)]
            for i in range(start, stop, increment):
                print(f'{algorithm}\t{i}', file = sys.stderr, end = '\t', flush = True)

                for _ in range(replicates):
                    print('r', end = '', file = sys.stderr, flush = True)    
                    
                    o = hmm_binding.binded_HMM(i, 4, hmmType = hmmType)
                    set_random(o)

                    t0 = time.time()
                    getattr(o, algorithm)(test_standard_data, **kwargs)
                    t1 = time.time()
                    o.deallocate()

                    print(f'statespace, {i}, {t1-t0}, {algorithm}, {o.hmmType}, {algorithm_version}')
                print('', file = sys.stderr, flush = True) # newline


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

        ## RSB ##        
        standard_test_statespace("viterbi", "RSB", inputsize, start, stop, increment, file)        
        standard_test_statespace("posteriorDecoding", "RSB", inputsize, start, stop, increment, file)
        standard_test_statespace("forward", "RSB", inputsize, start, stop, increment, file)
        standard_test_statespace("backward", "RSB", inputsize, start, stop, increment, file) 
        standard_test_statespace("baumWelch", "RSB", inputsize, start, stop, increment, file, '1', n_iterations = 1)

        


    ## Varying sparseness ##
    if 'sparse' in given_arguments:
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
        replicates = 5
        inputsize = 1500
        file = '../../test_framework/data/pantro3_X.fasta'

        """
        standard_test_sparseness("viterbi", "Conventional", inputsize, hidden_states, start, stop, increment, file)
        standard_test_sparseness("posteriorDecoding", "Conventional", inputsize, hidden_states, start, stop, increment, file)
        standard_test_sparseness("forward", "Conventional", inputsize, hidden_states, start, stop, increment, file)
        standard_test_sparseness("backward_time", "Conventional", inputsize, hidden_states, start, stop, increment, file)
        standard_test_sparseness("baumWelch", "Conventional", inputsize, hidden_states, start, stop, increment, file, str(1), n_iterations = 1)
        

        ## BLAS #
        standard_test_sparseness("viterbi", "BLAS", inputsize, hidden_states, start, stop, increment, file)
        standard_test_sparseness("posteriorDecoding", "BLAS", inputsize, hidden_states, start, stop, increment, file)
        standard_test_sparseness("forward", "BLAS", inputsize, hidden_states, start, stop, increment, file)
        standard_test_sparseness("backward_time", "BLAS", inputsize, hidden_states, start, stop, increment, file)
        standard_test_sparseness("baumWelch", "BLAS", inputsize, hidden_states, start, stop, increment, file, str(1), n_iterations = 1)
        

        ## CSR #
        standard_test_sparseness("viterbi", "CSR", inputsize, hidden_states, start, stop, increment, file)
        standard_test_sparseness("posteriorDecoding", "CSR", inputsize, hidden_states, start, stop, increment, file)
        standard_test_sparseness("forward", "CSR", inputsize, hidden_states, start, stop, increment, file)
        standard_test_sparseness("backward_time", "CSR", inputsize, hidden_states, start, stop, increment, file)
        standard_test_sparseness("baumWelch", "CSR", inputsize, hidden_states, start, stop, increment, file, str(1), n_iterations = 1)"""

        for hidden_states in [8]:
        
            print('# hs:', hidden_states, file = sys.stderr)
            ## Conventional #

            standard_test_sparseness("forward", "Conventional", inputsize, hidden_states, start, stop, increment, file, hidden_states)
            standard_test_sparseness("backward_time", "Conventional", inputsize, hidden_states, start, stop, increment, file, hidden_states)

            ## BLAS #
            standard_test_sparseness("forward", "BLAS", inputsize, hidden_states, start, stop, increment, file, hidden_states)
            standard_test_sparseness("backward_time", "BLAS", inputsize, hidden_states, start, stop, increment, file, hidden_states)

            ## CSR #
            standard_test_sparseness("forward", "CSR", inputsize, hidden_states, start, stop, increment, file, hidden_states)
            standard_test_sparseness("backward_time", "CSR", inputsize, hidden_states, start, stop, increment, file, hidden_states)

            ## RSB #
            standard_test_sparseness("forward", "RSB", inputsize, hidden_states, start, stop, increment, file, hidden_states)
            standard_test_sparseness("backward_time", "RSB", inputsize, hidden_states, start, stop, increment, file, hidden_states)

            
            


    ## Varying alphabet size ##
    if 'alphabet' in given_arguments:

        def standard_test_alphabet(algorithm, hmmType, input_size, start, stop, increment, file= '', algorithm_version = '', **kwargs):
            """ This standard test tests a varying size alphabet with a constant statespace and inputsize."""

            for i in range(start, stop, increment):
                print(f'{algorithm}\t{i}', file = sys.stderr, end = '\t', flush = True)
                test_standard_data = random.choices([j for j in range(i)], k = input_size*60) # Generates a data set with an arbitrary alphabet size (uniform).

                for _ in range(replicates):
                    print('r', end = '', file = sys.stderr, flush = True)
                    
                    o = hmm_binding.binded_HMM(7, i, hmmType = hmmType)
                    set_random(o)

                    t0 = time.time()
                    #test_standard_output = getattr(o, algorithm)(test_standard_data, **kwargs)
                    getattr(o, algorithm)(test_standard_data, **kwargs)
                    t1 = time.time()
                    o.deallocate()

                    print(f'alphabetsize, {i}, {t1-t0}, {algorithm}, {o.hmmType}, {algorithm_version}')
                print('', file = sys.stderr, flush = True)
                
        print('## Testing varying alphabet size ##', file = sys.stderr)
        inputsize = 1500 # the input size is constant. This number will be multiplied with 60 to become relatable with the other tests
        start = 2 # the state space
        stop = 15 
        increment = 2
        replicates = 5
        file = '../../test_framework/data/pantro3_X.fasta'

        ## Conventional ##
        standard_test_alphabet("viterbi", "Conventional", inputsize, start, stop, increment, file)
        standard_test_alphabet("posteriorDecoding", "Conventional", inputsize, start, stop, increment, file)
        standard_test_alphabet("forward", "Conventional", inputsize, start, stop, increment, file)
        standard_test_alphabet("backward", "Conventional", inputsize, start, stop, increment, file)
        standard_test_alphabet("baumWelch", "Conventional", inputsize, start, stop, increment, file, '1', n_iterations = 1)

        ## BLAS ##
        standard_test_alphabet("viterbi", "BLAS", inputsize, start, stop, increment, file)
        standard_test_alphabet("posteriorDecoding", "BLAS", inputsize, start, stop, increment, file)
        standard_test_alphabet("forward", "BLAS", inputsize, start, stop, increment, file)
        standard_test_alphabet("backward", "BLAS", inputsize, start, stop, increment, file)
        standard_test_alphabet("baumWelch", "BLAS", inputsize, start, stop, increment, file, '1', n_iterations = 1)

        ## CSR ##
        standard_test_alphabet("viterbi", "CSR", inputsize, start, stop, increment, file)
        standard_test_alphabet("posteriorDecoding", "CSR", inputsize, start, stop, increment, file)
        standard_test_alphabet("forward", "CSR", inputsize, start, stop, increment, file)
        standard_test_alphabet("backward", "CSR", inputsize, start, stop, increment, file)
        standard_test_alphabet("baumWelch", "CSR", inputsize, start, stop, increment, file, '1', n_iterations = 1)

        ## RSB ##        
        standard_test_alphabet("viterbi", "RSB", inputsize, start, stop, increment, file)        
        standard_test_alphabet("posteriorDecoding", "RSB", inputsize, start, stop, increment, file)
        standard_test_alphabet("forward", "RSB", inputsize, start, stop, increment, file)
        standard_test_alphabet("backward", "RSB", inputsize, start, stop, increment, file) 
        standard_test_alphabet("baumWelch", "RSB", inputsize, start, stop, increment, file, '1', n_iterations = 1)

        
        
