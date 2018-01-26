import numpy as np
import math

_P4VECTOR__constant_c = 3.0e8

class P4VECTOR():

    @staticmethod
    def AMATRIX(beta):
        if(abs(beta)>=1):
            print "************************"
            print "what the heck r u doing?"
            print "************************"
            raise Exception("I know python!")           
        A = [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]
        
        gamma = 1/math.sqrt(1-beta*beta)
        A[0][0] = A[3][3] = gamma
        A[1][1] = A[2][2] = 1
        A[0][3] = A[3][0] = -gamma*beta
        #print beta,gamma, A
        return np.matrix(A)


    def __init__(self,*args):
        
        if(not(len(args)==2  and args[1]<__constant_c) and len(args)!=0 and len(args)!=4):
            print "************************"
            print "what the heck r u doing?"
            print "************************"
            raise Exception("I know python!")
        
        if(len(args)==0):
            self.array = np.matrix([0,0,0,0]).getT()
            return
        
        if(len(args)==4):
            self.array = np.matrix(args).getT()
            return

        m = args[0]
        v0 = args[1]

        beta_tmp = v0/__constant_c
        gamma_tmp = 1/math.sqrt(1-beta_tmp*beta_tmp)
        Ep = m* __constant_c *gamma_tmp
        p3 = m*v0*gamma_tmp
        self.array = np.matrix([Ep,0,0,p3]).getT()

    def __str__(self):
        return str(self.array)

    def boost(self,A):
        self.array = np.dot(A,self.array)
        return self

#A = P4VECTOR.AMATRIX(0.1)
#p4 =  P4VECTOR(1,0)
#print p4.boost(A)
