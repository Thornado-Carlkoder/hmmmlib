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
    hmm_obj.setEmissionProbs(random_matrix(hmm_obj.n_hiddenstates, hmm_obj.n_observations, sparseness = 0))
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


