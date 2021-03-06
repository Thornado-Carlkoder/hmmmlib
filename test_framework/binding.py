import ctypes as c
import os

# authors: Thornado & Carl Koder

"""
Description:

This file creates a mirrored version of the HMM struct from hmm.c in the c-library
In then mirrors all its methods and algorithms in python such that all functionality in the c-library is accessible from python.

This file is used in test.py, which validates output of all the algorithms
It is also used by running_time.py, which calls each function on increasing parameters and records the timing.
"""

class HMM(c.Structure):
    """ creates a struct to match HMM """
    _fields_ = [("hiddenStates", c.c_uint),
                ("observations", c.c_uint),
                ("transitionProbs", c.POINTER(c.c_double)),
                ("emissionProbs", c.POINTER(c.c_double)),
                ("initProbs", c.POINTER(c.c_double))]

class binded_HMM:

    def __init__(self, n_hiddenstates, n_observations, hmmType = None):
        
        
        address_to_so = os.path.dirname(os.path.realpath(__file__)) + "/../hmmmlib/build/libHMMLIB.so"
        # Load the shared library into ctypes.

        self.libhmm = c.CDLL(os.path.abspath(address_to_so))
        
        # Set restypes for internal functions.
        #self.libhmm.HMMCreate.restype = c.POINTER(HMM)
        self.libhmm.HMMConventional.restype = c.POINTER(HMM)
        self.libhmm.HMMConventionalsparse.restype = c.POINTER(HMM)
        self.libhmm.HMMBLAS.restype = c.POINTER(HMM)
        self.libhmm.HMMCsr.restype = c.POINTER(HMM)
        
        
        
        
        self.libhmm.validateHMM.restype = c.c_bool
        self.libhmm.printHMM.restype = c.c_void_p
        self.libhmm.HMMDeallocate.restype = c.c_void_p

        # Set restypes for algorithms.
        #self.libhmm.forward.restype = c.POINTER(c.c_double)
        self.libhmm.F.restype = c.c_void_p
        self.libhmm.B.restype = c.c_void_p
        self.libhmm.viterbi.restype = c.POINTER(c.c_uint)
        self.libhmm.baumWelch.restype = c.c_void_p
        self.libhmm.posteriorDecoding.restype = c.POINTER(c.c_uint)

        self.n_hiddenstates = n_hiddenstates
        self.n_observations = n_observations

        self.hmmType = hmmType

        # Create HMM object
        if hmmType is None:
            print("Warning, the hmmType defaults to Conventional when none is given.")

        if hmmType == "Conventional" or hmmType is None:
            self.hmm = self.libhmm.HMMConventional(n_hiddenstates, n_observations)
            #print(" (A conventional hmm was created)")
        elif hmmType == "Consparse" or hmmType == "Conventional sparse":
            self.hmm = self.libhmm.HMMConventionalsparse(n_hiddenstates, n_observations)
            #print(" (A conventional sparse hmm was created)")
        elif hmmType == "BLAS":
            self.hmm = self.libhmm.HMMBLAS(n_hiddenstates, n_observations)
            #print(" (A BLAS hmm was created)")
        elif hmmType == "CSR": # Compressed Sparse Row
            self.hmm = self.libhmm.HMMCsr(n_hiddenstates, n_observations)
        else:
            raise ValueError("The hmmType argument given ({hmmType}) to binded_HMM() is invalid. \
                Please give any of None, 'Conventional', 'BLAS', 'CSR' or 'RSB'. ")



    def presentHMM(self):
        
        

        print('Presenting the HMM with the presentHMM()-function from the python-binding')
        print(' hiddenStates =', self.n_hiddenstates)
        print(' observations =', self.n_observations)
        
        print()
        formattedInitProbs = ["{:7.3f}".format(self.hmm[0].initProbs[i]) for i in range(self.n_hiddenstates)]
        print(' initProbs: [self.n_hiddenstates]\n ', ''.join(formattedInitProbs))
        print('\t(' + str(sum([self.hmm[0].initProbs[i] for i in range(self.n_hiddenstates)])) + ')')


        print()
        print(' transitionProbs: [self.n_hiddenstates][self.n_hiddenstates]', end = '\n  ')
        for row in range(self.n_hiddenstates):
            row_sum = 0
            for col in range(self.n_hiddenstates):
                value = self.hmm[0].transitionProbs[row*self.n_hiddenstates+col]
                print("{:7.3f}".format(value), end = ' ')
                row_sum += value
            print(end = '\t(' + str(row_sum) + ')\n  ')
        print()

        
        print(' emissionProbs: [self.n_hiddenstates][self.n_observations]', end = '\n  ') # [7][4] eller [self.n_hiddenstates][self.n_observations]
        for row in range(self.n_hiddenstates):
            row_sum = 0
            for col in range(self.n_observations):
                value = self.hmm[0].emissionProbs[row*self.n_observations+col]
                print("{:7.3f}".format(self.hmm[0].emissionProbs[row*self.n_observations+col]), end = ' ')
                row_sum += value
            print(end = '\t(' + str(row_sum) + ')\n  ')
        print()
        print(' The internal validation state is:', self.libhmm.validateHMM(self.hmm))

    def validate(self):
        return self.libhmm.validateHMM(self.hmm)

    ## Setters ##
    def setInitProbs(self, pi):
        if len(pi) != self.n_hiddenstates:
            raise Exception('Failed to set initProbs[]. initProbs[] should contain {a}'\
            'values but {b} were given.'.format(a = self.n_hiddenstates, b = len(pi)))
            #raise Exception('error1')
            
        self.hmm[0].initProbs = (c.c_double * self.n_hiddenstates)(*pi)


    def setTransitionProbs(self, new_trans_p):
        if len(new_trans_p) != self.n_hiddenstates:
            raise Exception('Failed to set transitionProbs[]. transitionProbs[] should contain {a}'\
            ' rows but {b} were given.'.format(a = self.n_hiddenstates, b = len(new_trans_p)))
            #raise Exception('error2')
            
        for row in new_trans_p:
            if len(row) != self.n_hiddenstates:
                raise Exception('Failed to set transitionProbs[]. transitionProbs[] should contain {a}'\
                ' columns but {b} were given.'.format(a = self.n_hiddenstates, b = len(row)))
                #raise Exception('error3')
                

        one_dimensional = [j for sub in new_trans_p for j in sub]
        # print(one_dimensional)
        self.hmm[0].transitionProbs = (c.c_double * (self.n_hiddenstates * self.n_hiddenstates))(*one_dimensional)


    def setEmissionProbs(self, new_emiss_p):
        if len(new_emiss_p) != self.n_hiddenstates:
            raise Exception('Failed to set emissionProbs[]. emissionProbs[] should contain {a}'\
            ' rows but {b} were given.'.format(a = self.n_hiddenstates, b = len(new_emiss_p)))
            #raise Exception('error4')
            
        for row in new_emiss_p:
            if len(row) != self.n_observations:
                raise Exception('Failed to set emissionProbs[]. emissionProbs[] should contain {a}'\
                ' columns but {b} were given.'.format(a = self.n_observations, b = len(row)))
                #raise Exception('error5')
                

        one_dimensional = [j for sub in new_emiss_p for j in sub]
        # print(one_dimensional)
        self.hmm[0].emissionProbs = (c.c_double * (self.n_hiddenstates * self.n_observations))(*one_dimensional)


    ## Getters ##
    def getInitProbs(self):
        return [self.hmm[0].initProbs[i] for i in range(self.n_hiddenstates)]

    def getTransitionProbs(self):
        self.n_hiddenstates = self.n_hiddenstates
        return [[self.hmm[0].transitionProbs[row * self.n_hiddenstates + col] for col in range(self.n_hiddenstates)] for row in range(self.n_hiddenstates)]

    def getEmissionProbs(self):
        self.n_hiddenstates = self.n_hiddenstates
        self.n_observations = self.n_observations
        return [[self.hmm[0].emissionProbs[row*self.n_observations + col] for col in range(self.n_observations)] for row in range(self.n_hiddenstates)]
        

    ## Algorithms ##
    def forward(self, observation_data):
        """ Returns a tuple. 1: pointer to alpha 2: Pointer to scalefactor """
        
        # Allocate scalefactor
        scalefactor = len(observation_data) * [0]
        scalefactor_c = (c.c_double * len(observation_data))(*scalefactor)

        # Allocate alpha matrix
        alpha_matrix = len(observation_data) * self.n_hiddenstates * [0]
        alpha_matrix_c = (len(observation_data) * self.n_hiddenstates * c.c_double)(*alpha_matrix)
        
    
        self.libhmm.F(self.hmm,
                      (c.c_int * len(observation_data))(*observation_data),
                      len(observation_data),
                      scalefactor_c,
                      alpha_matrix_c)
        
        return alpha_matrix_c, scalefactor_c
    
    def backward(self, observation_data, scalefactor = None):
        """ Returns a tuple. 1: pointer to alpha 2: Pointer to scalefactor
            time_test_only overrides scalefactor"""

        # Automatically calculate scalefactor with forward, if it is missing.
        if scalefactor is None: # Compute scalefactor yourself
            scalefactor = self.forward(observation_data)[1]
        
        # Allocate beta matrix
        beta_matrix = len(observation_data) * self.n_hiddenstates * [0]
        beta_matrix_c = (len(observation_data) * self.n_hiddenstates * c.c_double)(*beta_matrix)
        
    
        self.libhmm.B(self.hmm,
                      (c.c_int * len(observation_data))(*observation_data),
                      len(observation_data),
                      scalefactor,
                      beta_matrix_c)
        return beta_matrix_c, scalefactor # Returning the scalefactor might not be necessary, but makes it easier to handle and check output.

    
    def backward_time(self, observation_data):
        # This is more or less a copy of backward() with the slight modification that it doesn't care at all about the scale factor matrix.
        """ Returns a tuple. 1: pointer to alpha 2: Pointer to scalefactor
            time_test_only overrides scalefactor"""
    
        # Allocate unit scale factor (dummy for time tests)
        unit_vector = len(observation_data) * [1]
        scalefactor = (c.c_double * len(unit_vector))(*unit_vector)
        
        # Allocate beta matrix
        beta_matrix = len(observation_data) * self.n_hiddenstates * [0]
        beta_matrix_c = (len(observation_data) * self.n_hiddenstates * c.c_double)(*beta_matrix)
        
    
        self.libhmm.B(self.hmm,
                      (c.c_int * len(observation_data))(*observation_data),
                      len(observation_data),
                      scalefactor,
                      beta_matrix_c)
        return beta_matrix_c, scalefactor # Returning the scalefactor might not be necessary, but makes it easier to handle and check output.



    
    

    def viterbi(self, observation_data):
        output = (c.c_uint * len(observation_data))(*([0]*len(observation_data)))
        dummy_output = self.libhmm.viterbi(self.hmm,
                                     (c.c_uint * len(observation_data))(*observation_data),
                                     len(observation_data),
                                     output)
        return [output[i] for i in range(len(observation_data))] # Evt. generator?


    def posteriorDecoding(self, observation_data):
        output = (c.c_uint * len(observation_data))(*([0]*len(observation_data)))
        dummy_output = self.libhmm.posteriorDecoding(self.hmm,
                                               (c.c_uint * len(observation_data))(*observation_data),
                                               len(observation_data),
                                               output)
        return [output[i] for i in range(len(observation_data))] # Evt. generator?


    
    
    def baumWelch(self, observation_data, n_iterations = 1):
        _ = self.libhmm.baumWelch(self.hmm,
                                  (c.c_int * len(observation_data))(*observation_data),
                                  len(observation_data),
                                  n_iterations)
        return True


    def deallocate(self):
        c_struct = c.POINTER(HMM)(self.hmm)
        self.libhmm.HMMDeallocate(c_struct)

    def __del__(self):
        """ called when self’s reference count reaches zero. """
        self.deallocate()


