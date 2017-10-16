import re
from scipy.interpolate import interp1d
import numpy as np

class PowerTable(object):
    def __init__(self,fname):

        sci = '(\d+\.\d+E-*\d+)'
        flo = '(\d+\.\d+)'
        unit = '(\w+)'
        fmt = ('\s*'+flo+'\s'+unit+'\s+'+sci+'\s+'+sci
                +'\s+'+flo+'\s+'+unit
                +'\s+'+flo+'\s+'+unit
                +'\s+'+flo+'\s+'+unit
                +'.*')

        data = []

        ## MeV:1, cm:1
        units = {'keV':0.001,'MeV/mm':10,'MeV':1}

        with open(fname) as f:
            for line in f:
                res = re.match(fmt,line)
                if res is None: continue
                E,Eu, dE_dX1, dE_dX2, L = res.groups()[:5]
                data.append({'E':float(E)*units[Eu],'Range':float(L)/10.0})

	if not data:
            raise Exception("can you fucking at least get a valid table?")

        self.MaxR = max(_['E'] for _ in data)
        self.MaxE0 = max(_['Range'] for _ in data)

        self.R2E = interp1d(*zip(*[(_['Range'],_['E'],) for _ in data]),
                fill_value='extrapolate' )
        self.E2R = interp1d(*zip(*[(_['E'],_['Range'],) for _ in data]),
                fill_value='extrapolate' )

    def GetE(self,E0,X):
        if np.logical_or.reduce(np.array(E0)>self.MaxE0):
            raise Exception("exceed the maximum energy range of the table")

        R = self.E2R(E0)
        R_ = R-X
        E = self.R2E(R_)
        return np.clip(E,0,float('inf'))

    def GetE0(self,x):
        if np.logical_or.reduce(np.array(x)>self.MaxR):
            raise Exception("exceed the maximum distance range of the table")

        return self.R2E(x)
