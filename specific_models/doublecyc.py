"""
Determine the Tafel slope and pH dependence of a cycle with multiple parallel
pathways with different intermediate energies.
We saw the half order behavior when the intermediate energies were spread over a distribution,
so including additional parallel pathways should converge.
"""

from tafel import ElecMech
from tafel import SimTafel
import numpy as np
import sys

class MultiCycle(ElecMech.ElecMech):
    """
    Set a series of parallel cycles w/ individual rates and concentration dependencies
    """

    def __init__(self, cycle_tuple, weights=None, conc=None, dG=0., a=0.5):
        """
        Need to have initialized all ElecMech in cycle_tuple with rate constants beforehand
        -----------
        Parameters:
        cycle_tuple: (tuple) of ElecMech classes that act as parallel catalytic cycles
        weights: (list) population of each cycle (normalized to 1). fraction of active sites for each mech
        conc: pH, O2 concentration, depending on ElecMechs
        dG: (float) free energy of overall cycle
        -------
        Returns:
        None
        """

        super().__init__(dG=dG, a=a)
        self.cyc_tuple = cycle_tuple
        if weights is None:
            weight = 1.0/float(len(self.cyc_tuple)) #want the weights to add to 1 for comparison to indiv cycles
            self.weights = [weight]*len(self.cyc_tuple)
        else:
            self.weights = weights

    @property
    def pH(self):
        return self._pH

    @pH.setter
    def pH(self, value):
        for cyc in self.cyc_tuple:
            cyc._pH = value
        self._pH = value

    def rate(self, V):
        rate = 0
        for i, cyc in enumerate(self.cyc_tuple):
            rate += cyc.rate(V) * self.weights[i]
        return rate

def dG2k(dG_list):
    ##Generate rate constants from intermediate free energies
    kT = 8.61733*10**-5 * 298.15 #in eV/K, assuming T = 298.15
    ks = []
    for dG in dG_list:
        k = np.exp(-0.5*dG/kT)
        kn = np.exp(0.5*dG/kT)
        ks.append(k)
        ks.append(kn)
    return ks

def CycDistribution(em_cyc, delta_arr, dG_list, cyc_count=1, conc=[1.,1.]):
    cyc_list = []
    for i in range(cyc_count):
        dG1 = dG_list[0] + delta_arr[i]
        dG2 = dG_list[1] - delta_arr[i]
        ks = dG2k([dG1, dG2])
        cyc_list.append(em_cyc(k_rate=ks, conc=conc))
    return cyc_list

def genMCrates(mc, pH_range, V=-1.2):
    rate_list_mc = []
    for pH in pH_range:
        mc.pH = pH
        r = mc.rate(V)
        rate_list_mc.append(np.log10(np.abs(mc.rate(V))))
    return rate_list_mc

def plotMC(pH_range, rate_list_mc, dlogRdpH, ax0, ax1):
    ax0.plot(pH_range, rate_list_mc, linewidth=2)
    ax1.plot(pH_range, dlogRdpH, linewidth=2)

def plotParallelCycles(dGs, delta, pH_range, cyc_count, ax0, ax1, cyc_weights=None):

    if type(delta) == float:
        if cyc_count == 1:
            delta_arr = np.zeros(1)
        else:
            delta_arr = np.linspace(-0.5*delta, 0.5*delta, cyc_count)
    elif type(delta) == np.ndarray:
        delta_arr = delta

    cyc_list = CycDistribution(ElecMech.Rev2PcetO2, delta_arr, dGs, cyc_count=cyc_count)

    mc = MultiCycle(cyc_list, weights=cyc_weights)
    rates = genMCrates(mc, pH_range)
    dlogRdpH = np.gradient(np.array(rates), pH_range[1]-pH_range[0])
    plotMC(pH_range, rates, dlogRdpH, ax0, ax1)

if __name__ == "__main__":

    import matplotlib as mpl
    import matplotlib.pyplot as plt

    font = {'size':22}
    mpl.rc('font',**font)

    cyc_count = 1
    delta = 0.1
    sig = delta/2.
    V = -1.2
    pH_range = np.linspace(-15.,15.,400)
    dG = -0.1 # reaction is overall downhill

    dG1 = 0.2 #free energy barrier for reaction to intermediate 1
    dG2 = dG - dG1 # dG2 is defined implicitly by the overall energy of the cycle
    dGs = [dG1, dG2]

    #ks1 = dG2k([dG1, dG2])
    #cyc1 = ElecMech.Rev2PcetO2(k_rate=ks1, conc=[1.,1.])

    #sim = SimTafel.SimTafel(mc) #not currently using any methods from this class

    fig, (ax0,ax1) = plt.subplots(nrows=2)

    for i in range(1, cyc_count+1,8):
        #delta, weight_dist = np.polynomial.hermite.hermgauss(i)  #gaussian quadrature, x=deltas, y=weights
        #weight_dist *= 1./(np.pi**0.5) # the weights from gaussian quadrature actually sum to 2
        #delta *= sig*2**0.5
        weight_dist = np.ones(i)/i
        plotParallelCycles(dGs, delta, pH_range, i, ax0, ax1, cyc_weights=weight_dist)

    #ax1.set_ylim([-2,2])
    ax0.set_ylabel("Reaction Rate")
    ax1.set_ylabel("dR/dpH")
    plt.xlabel('pH')

    fig.set_size_inches(11.,11.,forward=True)

    plt.show()
    #plt.savefig("DoubleCycle_PcetO2_CycDist5_gaussDistD3.png", transparent=True, bbox_inches='tight', pad_inches=0.02)
