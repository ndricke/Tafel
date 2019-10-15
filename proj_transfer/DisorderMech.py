import sys

import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.integrate as integrate

import ElecMech
import SimTafel

font = {'size':22}
mpl.rc('font',**font)

"""
TODO:
1. Implement disorder in both rate constants and perceived voltage (and test effect of each)
2. Use gaussian quadrature for integrating over disorder
3. Write functions for fitting rate constants and disorder to experimental data
"""

class DisorderMech(ElecMech.ElecMech):
    """Model electrochemical cycles with disorder
    Note that there are several distinct ways of including disorder:
    1. Disorder in the voltage felt at different sites
    2. Disorder in the rate constants (need to ensure this is consistent for fwd and reverse)
    """

    def __init__(self, ks, sig):
        """
        Parameters:
        -----------
        ks: List
        rate constants for catalytic cycle
        sig: Float
        metric for disorder in the catalysis process

        Returns:
        --------
        None
        """
        ElecMech.ElecMech.__init__(self, k_rate=ks) #set general & rate constants
        self._sig = sig #energy variance of intermediate energies on surface
        self.V = None

    @property
    def sig(self):
        return self._sig

    @sig.setter
    def sig(self, value):
        self.sig2 = value*value
        dEspan = self.sig * 10 + self.sig2*200 #so the numerical integration works over a range of values of sig
        self.dE_range = [-1.*dEspan, dEspan]
        self._sig = value

    def probE(self, dE):
        """Energy probability distribution for surface sites assumed gaussian"""
        return np.exp(-dE**2/(2*self.sig2))/(2*np.pi*self.sig2)**0.5

    def disorderRate(self, dE):
        """Calculate rate for a given dE, assuming constant V"""
        return self.baseRate(self.V, dE)*self.probE(dE)

    def rate(self, V):
        """Calculate rate integrated over dE, the range of intrinsic energy disorder of surface"""
        self.V = V
        return integrate.quad(self.disorderRate, self.dE_range[0], self.dE_range[1])[0]

class DisorderPEET(DisorderMech):
    """Disorder included for PCET followed by electron transfer"""
    def __init__(self, ks, sig):
        super().__init__(ks, sig)

    def setConsts(self, ks):
        ## Set all of the base rate constants (excluding V, pH dependence)
        self.k1 = ks[0]
        self.kn1 = ks[1]
        self.k2 = ks[2]
        self.kn2 = ks[3]

    def baseRate(self, V, dE=0.):
        k1, kn1 = self.PCET(ks=[self.k1, self.kn1], V=V+dE, H=self.H, mech='acid')
        k2, kn2 = self.ET(ks=[self.k2, self.kn2], V=V-dE)
        return self.rev2(k1,k2,kn1,kn2)

class DisorderGCC(DisorderMech):
    """Disorder included for PCET followed by electron transfer"""
    def __init__(self, ks, sig):
        super().__init__(ks, sig)

    def setConsts(self, ks):
        ## Set all of the base rate constants (excluding V, pH dependence)
        self.k1 = ks[0]
        self.kn1 = ks[1]
        self.k2 = ks[2]
        self.kn2 = ks[3]

    def baseRate(self, V, dE=0.):
        k1, kn1 = self.PCET(ks=[self.k1, self.kn1], V=V+dE, H=self.H, mech='acid')
        k2, kn2 = self.PCET(ks=[self.k2, self.kn2], V=V-dE, H=self.O2, mech='acid')
        return self.rev2(k1,k2,kn1,kn2)



class DisorderHER(DisorderMech):

    def __init__(self, ks, sig):
        super().__init__(ks, sig)

    def setConsts(self, ks):
        ## Set all of the base rate constants (excluding V, pH dependence)
        self.k1 = ks[0]
        self.kn1 = ks[1]
        self.k2 = ks[2]
        self.kn2 = ks[3]

    def baseRate(self, V, dE=0.):
        k1, kn1 = self.PCET(ks=[self.k1, self.kn1], V=V+dE, H=self.H, mech='acid')
        k2, kn2 = self.PCET(ks=[self.k2, self.kn2], V=V-dE, H=self.H, mech='acid')
        return self.rev2(k1,k2,kn1,kn2)


if __name__ == "__main__":
    V_dom = np.linspace(-0.2,0.2,500) # domain for applied voltage
    dE_range = [-4,4]
    sig_list = [0.01, 0.05, 0.1, 0.2, 0.4] # list of values for varying disorders

    ## k1, kn1, k2, kn2
    k_list = [1,1000,1000,1]

    dis = DisorderHER(k_list, sig)
    rate_list = []
    disorder_list = []
    dis.H = 1.

    sim = SimTafel.SimTafel(dis)

    #dis.sig = 0.4
    #print(integrate.quad(dis.probE, dE_range[0], dE_range[1])) #check normalization of prob distribution
    #print(integrate.quad(dis.probE, dis.dE_range[0], dis.dE_range[1])) #check normalization of prob distribution
    #print(dis.rate(V=0.01))

    ## Plots current vs voltage for varying disorder
    for V in V_dom:
        rate_list.append(dis.rate(V))
    plt.scatter(V_dom, rate_list, label="0")

    for s in sig_list:
        dis.sig2 = s*s
        disorder_list = []
        for V in V_dom:
            disorder_list.append(dis.calcDisorderRate(V, dE_range)[0])
        plt.scatter(V_dom, disorder_list, label=str(s))

    print(len(rate_list))
    print(len(V_dom))

    plt.ylim([-5,5])
    plt.xlim([-0.1,0.1])
    plt.xlabel("Voltage (V)")
    plt.ylabel("Current")
    plt.legend()
    plt.show()



