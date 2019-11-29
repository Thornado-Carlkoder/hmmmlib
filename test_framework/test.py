from binding import *
from copy import copy

## These tests are equivalent to the ones in main_test.c
## Made in order to check that the python-binding is successfull.

"""
Description

This test file compares expected output with computed output for all algorithms. That means, that running 
this file successfully is proof that "everything works" and that the binding is sane.

This test file should be run from the hmmmlib/build/ directory:
    example:
    $ python ../../test_framework/test.py

"""



def test_viterbi(hmmType = "Conventional"):
    # create HMM for testing viterbi
    o = binded_HMM(7, 4, hmmType = hmmType)
    o.setInitProbs([0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00])
    o.setTransitionProbs([[0.00, 0.00, 0.90, 0.10, 0.00, 0.00, 0.00],
                        [1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
                        [0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00],
                        [0.00, 0.00, 0.05, 0.90, 0.05, 0.00, 0.00],
                        [0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00],
                        [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00],
                        [0.00, 0.00, 0.00, 0.10, 0.90, 0.00, 0.00]])
    o.setEmissionProbs([[0.30, 0.25, 0.25, 0.20],
                        [0.20, 0.35, 0.15, 0.30],
                        [0.40, 0.15, 0.20, 0.25],
                        [0.25, 0.25, 0.25, 0.25],
                        [0.20, 0.40, 0.30, 0.10],
                        [0.30, 0.20, 0.30, 0.20],
                        [0.15, 0.30, 0.20, 0.35]])

    test_viterbi_data = [3, 0, 2, 0, 2, 3, 3, 2, 3, 3, 2, 1, 3, 1, 0, 2, 3, 1, 0, 1, 0, 1, 0, 1, 3, 0, 1, 0, 2, 1, 0, 3, 0, 0, 0, 2, 3, 2, 0, 3, 0, 2, 3, 0, 2, 3, 0, 2, 2, 0, 1, 0, 2, 1, 0, 0, 3, 1, 2, 3, 1, 0, 3, 1, 0, 3, 1, 2, 0, 1, 0, 2, 0, 2, 1, 3, 0, 1, 2, 0, 3, 1, 0, 2, 3, 0, 2, 0, 0, 2, 0, 1, 0, 3, 1, 0, 3, 0, 0]
    test_viterbi_output = o.viterbi(test_viterbi_data)
    test_viterbi_expected = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1]
    for i, j in zip(test_viterbi_output, test_viterbi_expected):
        assert i == j
    
    returnvalue = copy(o.validate())
    o.deallocate()
    return returnvalue



def test_posterior_decoding(hmmType = None):
    # create HMM for testing posterior_decoding
    o = binded_HMM(7, 4, hmmType = hmmType)
    o.setInitProbs([0.00, 0.00, 0.00, 1.00, 0.00, 0.00, 0.00])
    o.setTransitionProbs([[0.00, 0.00, 0.90, 0.10, 0.00, 0.00, 0.00],
                        [1.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00],
                        [0.00, 1.00, 0.00, 0.00, 0.00, 0.00, 0.00],
                        [0.00, 0.00, 0.05, 0.90, 0.05, 0.00, 0.00],
                        [0.00, 0.00, 0.00, 0.00, 0.00, 1.00, 0.00],
                        [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 1.00],
                        [0.00, 0.00, 0.00, 0.10, 0.90, 0.00, 0.00]])
    o.setEmissionProbs([[0.30, 0.25, 0.25, 0.20],
                        [0.20, 0.35, 0.15, 0.30],
                        [0.40, 0.15, 0.20, 0.25],
                        [0.25, 0.25, 0.25, 0.25],
                        [0.20, 0.40, 0.30, 0.10],
                        [0.30, 0.20, 0.30, 0.20],
                        [0.15, 0.30, 0.20, 0.35]])

    test_posterior_decoding_data = [0, 1, 3, 3, 2, 1, 2, 3, 0, 0, 1, 2, 1, 0, 2, 3, 0, 1, 3, 1, 2, 0, 3, 1, 2, 0, 3, 1, 2, 3, 0, 1, 2, 1, 0, 3, 2, 1, 0, 0, 1, 0, 3, 2, 0, 3, 0, 0, 0, 3, 1, 1, 0, 0, 3, 1, 0, 3, 1, 0, 3, 0, 1, 0, 3, 2, 0, 1, 0, 3, 1, 0, 2, 1, 3, 0, 2, 1, 0, 3, 2, 0, 1, 2, 0, 3, 1, 0, 2, 3, 0, 1, 2, 0, 1, 3, 0, 1, 2, 0, 3, 1, 2, 3, 0]
    test_posterior_decoding_output = o.posteriorDecoding(test_posterior_decoding_data)
    test_posterior_decoding_expected = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0, 2, 1, 0]
    for i, j in zip(test_posterior_decoding_output, test_posterior_decoding_expected):
        assert i == j

    returnvalue = copy(o.validate())
    o.deallocate()
    return returnvalue





def test_forward(hmmType = None):
    o2 = binded_HMM(2, 2, hmmType = hmmType)
    o2.setInitProbs([.2, .8])
    o2.setTransitionProbs([[0.5, 0.5],
                        [0.3, 0.7]])
    o2.setEmissionProbs([[.3, .7],
                        [.8, .2]])
    data = [0, 0, 0, 0, 0, 1, 1, 0, 0, 0]

    output, scalefactor_from_forward = o2.forward(data)
    output = [output[i] for i in range(len(data)*len(o2.getTransitionProbs()))]
    expected = [0.08571428571428569, 0.9142857142857143, 0.14832962138084632, 0.8516703786191537, 0.1557067218345664, 0.8442932781654335, 0.15658545239699356, 0.8434145476030065, 0.15669025947764703, 0.843309740522353, 0.6342802366859542, 0.3657197633140457, 0.7227357992924823, 0.2772642007075177, 0.23084308303092038, 0.7691569169690796, 0.1656531912620611, 0.834346808737939, 0.1577734760873247, 0.8422265239126754]
    for i, j in zip(output, expected):
        assert (i - j) < 0.00001
        #print(i, j)


    returnvalue = copy(o2.validate())
    o2.deallocate()
    return returnvalue


def test_backward(hmmType = None):
    o2 = binded_HMM(2, 2, hmmType = hmmType)
    o2.setInitProbs([.2, .8])
    o2.setTransitionProbs([[0.5, 0.5],
                        [0.3, 0.7]])
    o2.setEmissionProbs([[.3, .7],
                        [.8, .2]])
    data = [0, 0, 0, 0, 0, 1, 1, 0, 0, 0]

    output, scalefactor_from_backward = o2.backward(data)
    output = [output[i] for i in range(len(data)*2)]
    
    # oexpected is the expected values of the same hmm, but with wide alpha and beta tables instead of narrow
    #oexpected = [0.838486, 0.848495, 0.854859, 0.898964, 1.267584, 1.076550, 0.944879, 0.862481, 0.868282, 1.000000, 1.015142, 1.026387, 1.026767, 1.018758, 0.950282, 0.867237, 1.143683, 1.041273, 1.026152, 1.000000]
    expected = [0.838486, 1.015142, 0.848495, 1.026387, 0.854859, 1.026767, 0.898964, 1.018758, 1.267584, 0.950282, 1.076550, 0.867237, 0.944879, 1.143683, 0.862481, 1.041273, 0.868282, 1.026152, 1.000000, 1.000000]
    

    for i, j in zip(output, expected):
        assert (i - j) < 0.00001
        #print(i, j)


    returnvalue = copy(o2.validate())
    o2.deallocate()
    return returnvalue








def test_baumwelch(hmmType = None):
    o2 = binded_HMM(2, 2, hmmType = hmmType)
    o2.setInitProbs([.2, .8])
    o2.setTransitionProbs([[0.5, 0.5],
                        [0.3, 0.7]])
    o2.setEmissionProbs([[.3, .7],
                        [.8, .2]])
    test_for_back_bw_data = [0, 0, 0, 0, 0, 1, 1, 0, 0, 0]


    test_baumwelch_output = o2.baumWelch(test_for_back_bw_data, 1)

    epsilon = 0.000001
    assert o2.hmm[0].transitionProbs[0] - 0.43921478415604948 < epsilon
    assert o2.hmm[0].transitionProbs[1*o2.hmm[0].hiddenStates + 1] - 0.78554317773503979 < epsilon
    assert o2.hmm[0].emissionProbs[0] - 0.46160107308583781 < epsilon
    assert o2.hmm[0].emissionProbs[1*o2.hmm[0].observations + 1] - 0.084984433203479412 < epsilon


    #o2.presentHMM()
    returnvalue = copy(o2.validate())
    o2.deallocate()
    return returnvalue


print('Conventional forward', test_forward("Conventional"))
print('Conventional backward', test_backward("Conventional"))
print('Conventional baum-welch', test_baumwelch("Conventional")) 
print('Conventional posterior decoding', test_posterior_decoding("Conventional"))
print('Conventional viterbi', test_viterbi("Conventional"))

print('Conventional sparse fw', test_forward("Consparse"))
print('Conventional sparse bw', test_backward("Consparse"))
print('Conventional sparse baum-welch', test_baumwelch("Consparse")) 
print('Conventional sparse posterior decoding', test_posterior_decoding("Consparse"))
print('Conventional sparse viterbi', test_viterbi("Consparse"))

print('BLAS fw', test_forward("BLAS"))
print('BLAS bw', test_backward("BLAS"))
print('BLAS baum-welch', test_baumwelch("BLAS")) 
print('BLAS posterior decoding', test_posterior_decoding("BLAS"))
print('BLAS viterbi', test_viterbi("BLAS"))

print('CSR forward', test_forward("CSR"))
print('CSR backward', test_backward("CSR"))
print('CSR baum-welch', test_baumwelch("CSR")) 
print('CSR posterior decoding', test_posterior_decoding("CSR"))
print('CSR viterbi', test_viterbi("CSR"))

print('RSB forward', test_forward("RSB"))
print('RSB backward', test_backward("RSB"))
print('RSB baum-welch', test_baumwelch("RSB")) 
print('RSB posterior decoding', test_posterior_decoding("RSB"))
print('RSB viterbi', test_viterbi("RSB"))
