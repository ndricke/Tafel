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

class DisorderMech(ElecMech.ElecMech):

    def __init__(self, ks, sig):
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
    V_dom = np.linspace(-0.2,0.2,500)
    ## k1, kn1, k2, kn2
    k_list = [1,1000,1000,1] 
    sig = 0.01
    dE_range = [-4,4]

    dis = DisorderHER(k_list, sig)
    rate_list = []
    disorder_list = []
    dis.H = 1.
#    sig_list = [0.001, 0.01, 0.05, 0.1, 0.2, 0.4, 0.8]
    sig_list = [0.01, 0.05, 0.1, 0.2, 0.4]

    sim = SimTafel.SimTafel(dis)

    #dis.sig = 0.4
    #print(integrate.quad(dis.probE, dE_range[0], dE_range[1])) #check normalization of prob distribution
    #print(integrate.quad(dis.probE, dis.dE_range[0], dis.dE_range[1])) #check normalization of prob distribution
    #print(dis.rate(V=0.01))

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

    ## pH dependence with intrinsic disorder
    #"""
    V_list = []
    onset_J = 0.1
    n = 40
    pH_list = np.linspace(-2, 8, n)

    fig, (ax0,ax1) = plt.subplots(nrows=2)

    for sig in sig_list:
        dis.sig = sig
        V_arr = np.zeros(n)
        for i, pH in enumerate(pH_list):
            dis.pH = pH
            V_arr[i] = sim.findOnsetV(onset_J, onsetV_guess=-0.03)['x']

        ax0.plot(pH_list, V_arr)

        dVdpH = np.gradient(V_arr, pH_list[1]-pH_list[0]) #inverse of Tafel slope; calc this way bc V_dom indep var w/ even spacing
        ax1.plot(pH_list, dVdpH, label=str(sig))

    plt.xlabel("pH")
    ax0.set_ylabel("Onset Potential (V)")
    ax1.set_ylabel("dV/dpH")
    #plt.xlim([-0.01, 0.21])
    fig.set_size_inches(11.,11.,forward=True)
    plt.legend()
    #plt.savefig("IntrinsicHet_1t1000_pHvsOnsetV.png", transparent=True, bbox_inches='tight', pad_inches=0.02)
    plt.show()



    #"""

    ## Plot onset potential as a function of intrinsic disorder
    """
    V_list = []
    onset_J = 0.1
    n = 60
    sig_list = np.linspace(0.001, 0.2, n)

    fig, ax = plt.subplots()

    ## Test to make sure that disorder converges to ordered at sig = 0.
    #ordered_mech = ElecMech.RevPcet2()
    #ordered_mech.setConsts(k_list)
    #ordered_mech.setConc(0)
    #ordered_sim = SimTafel.SimTafel(ordered_mech)
    #V_ordered = ordered_sim.findOnsetV(onset_J)['x']
    #plt.scatter(0., V_ordered)

    k_ratios = [1, 5, 10, 25, 50, 100]
    for ratio in k_ratios:
        ks = [1./ratio,ratio,ratio,1./ratio]
        dis.setConsts(ks)
        print(dis.k1, dis.kn1, dis.k2, dis.kn2)
        V_arr = np.zeros(n)
        for i, sig in enumerate(sig_list):
            dis.sig = sig
            V_arr[i] = sim.findOnsetV(onset_J, onsetV_guess=-0.03)['x']

        plt.scatter(sig_list, V_arr, )
        plt.plot(sig_list, V_arr, label=ratio)

    plt.xlabel("Disorder")
    plt.ylabel("Onset Potential (V)")
    plt.xlim([-0.01, 0.21])
    fig.set_size_inches(11.,11.,forward=True)
    plt.legend()
    plt.savefig("IntrinsicHet_DisorderVsOnsetV.png", transparent=True, bbox_inches='tight', pad_inches=0.02)
    #plt.show()
    """


    ## Plot Tafel slope as a function of intrinsic disorder
    """

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
    fig.set_size_inches(11.,11., forward=True)

    #plt.show()
    plt.savefig("IntrinsicHet_1t1_VvsJ.png", transparent=True, bbox_inches='tight', pad_inches=0.02)
    """
