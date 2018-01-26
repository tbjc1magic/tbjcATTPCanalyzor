import math
import random
from P4VECTOR import P4VECTOR
_KINEMATICS__constant_u = 1.66e-27
_KINEMATICS__constant_MeV = 1.6e-13
_KINEMATICS__constant_e = 1.6e-19
_KINEMATICS__constant_c = 3.0e8

class KINEMATICS:

    def __init__(self,m, K0, Eex2, Eex3):
        """m:u,K0:MeV,Eex2:MeV,Eex3:MeV"""
        #print m,K0,Eex2,Eex3
        #######################
        ####input parameters###
        #######################

        self.K0 = K0
        self.m = m
        self.Eex2 = Eex2
        self.Eex3 = Eex3
        #######################
        ###output parameters###
        #######################

        self.K2 = None
        self.K3 = None
        self.thetalab2 = None
        self.philab2 = None
        self.thetalab3 = None
        self.philab3 = None
        self.V2 = None
        self.V3 = None
        self.thetaCMS = None
        self.phiCMS = None

    def reset(self):

        self.K2 = None
        self.K3 = None
        self.thetalab2 = None
        self.philab2 = None
        self.thetalab3 = None
        self.philab3 = None
        self.V2 = None
        self.V3 = None
        self.thetaCMS = None
        self.phiCMS = None

    def randomgenerate(self):
        cos_theta = random.random()*2-1
        phiCMS = random.random()*math.pi*2
        thetaCMS = math.acos(cos_theta)

        self.calculate(thetaCMS, phiCMS)

    def calculate(self, thetaCMS, phiCMS):
        """unit:thetaCMS(radians),phiCMS(radians)"""

        self.reset()

        m0_tmp = self.m[0]*__constant_u
        m1_tmp = self.m[1]*__constant_u
        m2_tmp = self.m[2]*__constant_u
        m3_tmp = self.m[3]*__constant_u

        Ek_tmp = self.K0*__constant_MeV
        v0_tmp = self.Ek2v(m0_tmp,Ek_tmp)
        P4Va = P4VECTOR(m0_tmp, v0_tmp)
        u_tmp = self.findCMSu(m0_tmp,m1_tmp,v0_tmp)
        beta_tmp = u_tmp/__constant_c
        A_tmp = P4VECTOR.AMATRIX(beta_tmp)
        _A_tmp = P4VECTOR.AMATRIX(-beta_tmp)
        P4Va.boost(A_tmp)

        P4Vb = P4VECTOR(m1_tmp,0)
        P4Vb.boost(A_tmp)
        #print P4Vb
        energyA_tmp = P4Va.array[0,0]*__constant_c
        energyB_tmp = P4Vb.array[0,0]*__constant_c

        DeltaE_tmp = (self.Eex2+self.Eex3)*__constant_MeV
        Etot_tmp = energyA_tmp + energyB_tmp - DeltaE_tmp

        E2_0_tmp = m2_tmp*__constant_c*__constant_c
        E3_0_tmp = m3_tmp*__constant_c*__constant_c

        tmp1 = E3_0_tmp*E3_0_tmp+E2_0_tmp*E2_0_tmp-Etot_tmp*Etot_tmp

        if(tmp1>0):
            print "************************"
            print "what the heck r u doing?"
            print "************************"
            raise Exception("I know shit!")

        tmp2 = 2*E3_0_tmp*E2_0_tmp

        Pcms2_tmp = (tmp1*tmp1-tmp2*tmp2)/(4*Etot_tmp*Etot_tmp)
        Pcms_tmp = math.sqrt(Pcms2_tmp)/__constant_c

        E3_tmp = math.sqrt(E3_0_tmp*E3_0_tmp+
                Pcms_tmp*Pcms_tmp*__constant_c*__constant_c)

        E2_tmp = math.sqrt(E2_0_tmp*E2_0_tmp+
                Pcms_tmp*Pcms_tmp*__constant_c*__constant_c)

        PcmsZ_tmp = Pcms_tmp*math.cos(thetaCMS)
        PcmsX_tmp = Pcms_tmp*math.sin(thetaCMS)*math.cos(phiCMS)
        PcmsY_tmp = Pcms_tmp*math.sin(thetaCMS)*math.sin(phiCMS)

        P4Vc = P4VECTOR(E2_tmp/__constant_c,PcmsX_tmp,PcmsY_tmp,PcmsZ_tmp)
        P4Vd = P4VECTOR(E3_tmp/__constant_c,-PcmsX_tmp,-PcmsY_tmp,-PcmsZ_tmp)

        P4Vc.boost(_A_tmp)
        P4Vd.boost(_A_tmp)

        K2_tmp = P4Vc.array[0,0]*__constant_c-m2_tmp*__constant_c*__constant_c
        K3_tmp = P4Vd.array[0,0]*__constant_c-m3_tmp*__constant_c*__constant_c

        K2 = K2_tmp/__constant_MeV
        K3 = K3_tmp/__constant_MeV

        P2_tmp = math.sqrt(math.pow(P4Vc.array[1,0],2)+
                            math.pow(P4Vc.array[2,0],2)+
                            math.pow(P4Vc.array[3,0],2))

        P3_tmp = math.sqrt(math.pow(P4Vd.array[1,0],2)+
                            math.pow(P4Vd.array[2,0],2)+
                            math.pow(P4Vd.array[3,0],2))

        thetalab2 = math.acos(P4Vc.array[3,0]/P2_tmp)
        philab2 = phiCMS

        thetalab3 = math.acos(P4Vd.array[3,0]/P3_tmp)
        if(phiCMS<math.pi):
            philab3 = phiCMS + math.pi
        else:
            philab3 = phiCMS - math.pi

        V2 = self.Ek2v(m2_tmp,K2_tmp)
        V3 = self.Ek2v(m3_tmp,K3_tmp)

        self.K2 = K2
        self.K3 = K3
        self.V2 = V2
        self.V3 = V3
        self.thetalab2 = thetalab2/math.pi*180
        self.philab2 = philab2/math.pi*180
        self.thetalab3 = thetalab3/math.pi*180
        self.philab3 = philab3/math.pi*180
        self.thetaCMS = thetaCMS/math.pi*180
        self.phiCMS = phiCMS/math.pi*180

        #print self.K2,self.K3, self.thetalab2, self.thetalab3

    def findCMSu(self, m1, m2, v0):
        """unit:m1(kg),m2(kg),v0(m/s)"""
        beta_tmp = v0 /__constant_c
        gamma_tmp = 1/math.sqrt(1-beta_tmp*beta_tmp)
        m1r = m1*gamma_tmp
        u = m1r*v0/(m1r+m2)
        return u

    def Ek2v(self, m, Ek):
        """unit:m(kg),v(m/s)"""
        ratio_tmp = Ek/(m*__constant_c*__constant_c)
        v=__constant_c*math.sqrt(1-1/((ratio_tmp+1)*(ratio_tmp+1)))
        return v
