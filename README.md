[![Build Status](https://travis-ci.com/Thornado-Carlkoder/hmmmlib.svg?branch=master)](https://travis-ci.com/Thornado-Carlkoder/hmmmlib)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/379148760d544cab8b4a14322400a1ea)](https://www.codacy.com/gh/Thornado-Carlkoder/hmmmlib?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=Thornado-Carlkoder/hmmmlib&amp;utm_campaign=Badge_Grade)



# HMMMLIB <img src='ressourcer/hmm_smaller.png' align="right" height="138.5" />

#### *Hidden Markov-Model Matrix Library*

* HMM library written in C 
* Has Python bindings allowing it to be used by Python as well
* The library should be made in a modular way such that the user can easily add functionalities


## GOAL

The goal of this project is:

	- To get a strong theoretical understanding of HMM's and different versions of the latter
	- Create a robust and easy to use library for both C and Python
	- When the project is finished, make it open source for everyone to add and use

### Project description:
The goal of the project is to implement an efficient C library for hidden Markov models that supports the Forward and Backward algorithms, Viterbi and Posterior decoding and parameter optimisation (Baum-Welch). The focus will be on algorithmic engineering for dense matrix HMMs and in particular we will compare the straightforward implementations with BLAS based linear algebra formulations of the algorithms.


### Learning goals:
* The student should be able to describe a general HMM, implement a HMM, the forward, backward, viterbi and posterior decoding algorithms - The student should be able to analyse the different implementations in terms of running time.
* The student should be able to discuss and evaluate the findings of the experiments in contrast to the decisions behind the implementation.


## WEEKLY MEETING:
* INTERNAL: Torsdag kl 11
* THOMAS: Torsdag kl 12

## GENERAL WORK PLAN
1. Make HMM with basic algorithms
2. Make a test framework with BLAS
3. Make sparse matrix optimized algorithms


### Evaluation problem
The evaluation problem regards the probability that a particular sequence of symbols is produced by a particular model.
For evaluation we use two algorithms: the forward algorithm or the backward algorithm

### Decoding problem
The decoding problem regards to inferring the most likely sequence of states produced by a given sequence of observables.
For decoding we apply the Viterbi and Posterior Decoding algorithms.

### Training problem
Training problem answers the question: Given a model structure and a set of sequences, find the model that best fits the data.
For this problem we can use the Baum Welch = forward\*backward algorithm

Possible optimisations: 

* Solving the "training" problem as a constrained optimisation problem and use lagrange multipliers?

## DOCUMENTATION

In this section all the documentation for the library is written as well as the different papers that have been implemented.

## TASK MANAGEMENT 

See the [GitHub project manager](https://github.com/Thornado-Carlkoder/hmmmlib/projects)



## RUNNING TIME TESTS

In order to quickly assess improvements to our algorithms, we have automated the running time tests.

The automated running time tests are part of the test framework - Thus they are located in the  test_framework/ directory.

Please note that density = 1 - sparseness 

These are the automated running time tests supported:

* Varying alphabet size: `running_time_alphabet_size.py`
* Varying input size: `running_time_input_size.py`
* Varying density: `running_time_sparseness.py`
* Varying state space: `running_time_state_space.py`
* Co-varying state space and density: `running_time_state_spaceVSsparseness.py`



Running each of these files with the python interpreter will run a series of tests. For instance, if we want to gauge the running time of increased input size, we will run the following command:

```sh
cd test_framework/
python running_time_input_size.py > input_size.csv
```

The scripts all output a csv-formatted file to STDOUT and progress diagnostics to STDERR. Because of this, it is recommended to pipe STDOUT to a separate file.



Note: If you're using Mac or Windows, the shared object file will likely have a varying file extension than what is hard coded in the binding.py file. Please change the `address_to_so` variable in binding.py.



If you want to change the paremeters passed to the algorithms, you should open the files individually and edit the hard coded variables. These variables are the following:

* `inputsize`: The length of the sequence that the HMM algorithm shold be applied to
* `start`, `stop` and `increment`: The parameter vector.
* `replicates`: The number of dependent replicates for each test. Used to infer the standard error of a measurement.
* `file`: The fasta-formatted file used as sequence in all running time tests but the alphabet size test.



## READING MATERIAL

* https://web.stanford.edu/class/cs262/archives/notes/lecture6.pdf
* https://en.wikipedia.org/wiki/Hidden_Markov_model

#### For c programming

* http://c-faq.com/aryptr/dynmuldimary.html
* https://www.dipmat.univpm.it/~demeio/public/the_c_programming_language_2.pdf
* https://developer.ibm.com/articles/au-toughgame/
* https://cdecl.org/
* http://derekmolloy.ie/hello-world-introductions-to-cmake/

### CMAKE

* http://derekmolloy.ie/hello-world-introductions-to-cmake/
