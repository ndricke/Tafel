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

class DisorderHER(ElecMech.ElecMech):

    def __init__(self, ks, sig):
        ElecMech.ElecMech.__init__(self, k_rate=ks) #set general & rate constants
        self._sig = sig #energy variance of intermediate energies on surface
        self.V = None

    def setConsts(self, ks):
        ## Set all of the base rate constants (excluding V, pH dependence)
        self.k1 = ks[0]
        self.kn1 = ks[1]
        self.k2 = ks[2]
        self.kn2 = ks[3]

    @property
    def sig(self):
        return self._sig

    @sig.setter
    def sig(self, value):
        self.sig2 = value*value
        dEspan = self.sig2 * 10
        self.dE_range = [-1.*dEspan, dEspan]
        self._sig = value


    def baseRate(self, V):
        k1, kn1 = self.PCET(ks=[self.k1, self.kn1], V=V, H=self.H, mech='acid')
        k2, kn2 = self.PCET(ks=[self.k2, self.kn2], V=V, H=self.H, mech='acid')
        return self.rev2(k1,k2,kn1,kn2)

    def probE(self, dE):
        """Energy probability distribution for surface sites assumed gaussian"""
        return np.exp(-dE**2/(2*self.sig2))/(2*np.pi*self.sig2)**0.5

    def disorderRate(self, dE):
        """Calculate rate for a given dE, assuming constant V"""
        return self.baseRate(self.V+dE)*self.probE(dE)

    def rate(self, V):
        """Calculate rate integrated over dE, the range of intrinsic energy disorder of surface"""
        self.V = V
        return integrate.quad(self.disorderRate, self.dE_range[0], self.dE_range[1])[0]


if __name__ == "__main__":
    V_dom = np.linspace(-0.2,0.2,500)
    k_list = [1,1,1,1] #k1, kn1, k2, kn2
    dE_range = [-10,10]
    sig = 0.01

    dis = DisorderHER(k_list, sig)
    rate_list = []
    disorder_list = []
    dis.H = 1.
    sig_list = [0.01, 0.05, 0.1, 0.2]

    ## Creates /work/Tafel/disorder/disorderHER_VvsI.png
    """
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
    """

    ##print(integrate.quad(dis.probE, dE_range[0], dE_range[1])) #check normalization of prob distribution

    sim = SimTafel.SimTafel(dis)

    print(dis.sig)
    dis.sig = 1
    print(dis.sig2)
    print(dis.sig)


    fig, (ax0,ax1) = plt.subplots(nrows=2)
    for s in sig_list:
        dis.sig = s
        logI_list = []; dV_list = [];
        for V in V_dom:
            I = dis.rate(V)
            logI = np.log10(np.abs(I))
            logI_list.append(logI)

        dlogIdV = np.gradient(logI_list, V_dom[1]-V_dom[0]) #inverse of Tafel slope; calc this way bc V_dom indep var w/ even spacing
        dVdlogI = 1./np.array(dlogIdV)

        ax0.plot(V_dom, logI_list)
        ax1.plot(V_dom, dVdlogI, label=str(s))

    plt.xlabel('V')
    #ax1.set_ylim(input_ylim)
    ax0.set_xlim([V_dom[0], V_dom[-1]])
    ax1.set_xlim([V_dom[0], V_dom[-1]])
    ax0.set_ylabel('log(I)')
    ax1.set_ylabel('$\partial$V/$\partial$log(I)')
    ax1.legend()


    plt.show()
