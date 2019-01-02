"""
Determine the Tafel slope and pH dependence of a cycle with 2 parallel pathways with different intermediate energies
We saw the half order behavior when the intermediate energies were spread over a distribution,
so including additional parallel pathways should converge


Tafel slope is of secondary interest to pH dependence, so we'll focus on the latter here first

The 2 requirements of SimTafel are a rate function that only requires the voltage, V, as input,
setConsts, which sets all the rate constants (k), and s.pH set as the same for all individual cycles

One obvious way to do this is to have a tuple of ElecMech objects, but currently each cycle has its own individual
concentration and voltage. I could potentially write routines to pass these to each of the cycles
--> setter for s.pH that also sets s.pH for each individual cycle
--> the rate function will call each cycle in order, finding individual rates, and return their sum
"""

from tafel import ElecMech
from tafel import SimTafel

class MultiCycle(ElecMech.ElecMech):

    def __init__(self, cycle_tuple, k_rate=None, conc=None, dG=0., a=0.5):
        super().__init__(k_rate, conc, dG, a)
        self.cyc_tuple = cycle_tuple

        ## Set a series of parallel cycles w/ individual rates and concentration dependencies

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
        for cyc in self.cyc_tuple:
            rate += cyc.rate(V)
        return rate

if __name__ == "__main__":

    import numpy as np
    import matplotlib.pyplot as plt


    dG = -0.1 # reaction is overall downhill 0.2

    kT = 8.61733*10**-5 * 298.15 #in eV/K, assuming T = 298.15

    k_ratio = np.exp(-1.*dG/(kT))

    #how to split k_ratio into k1/kn1 and k2/kn2? --> essentially need to assume barrier heights (possibly 0?)
    #if we assume the barrier height is 0
    #now, instead of using some ks2_speedup factor, maybe it would be better to just list 2 intermediate free energies
    dG1 = 0.2
    dG2 = 0.3

    #we want both cycles to have intermediates with the same barriers, changing only the intermediate energy
    k1_ratio = np.exp(-1.*dG1/(kT))
    k1 = 
    k2 =
    kn1 =
    kn2 =

    k2_ratio = np.exp(-1.*dG2/(kT)) # k1k2/knkn2 ratio for cycle 2

    #I think the barrier is determined by


    ks1 = [0.1,10.,10.,0.1] # k1, k2, kn1, kn2
    #ks2 = [0.1,10.,10.,0.1] # k1, k2, kn1, kn2
    ks2 = [1./11.,11.,11.,1./11.] # k1, k2, kn1, kn2

    #ks1 = [1.,1.,1.,1.] # k1, k2, kn1, kn2
    #ks2 = [1.,1.,1.,1.] # k1, k2, kn1, kn2

    pH_range = np.linspace(-7.,7.,400)
    V = 0.

    cyc1 = ElecMech.Rev2PcetO2(k_rate=ks1, conc=[1.,1.])
    cyc2 = ElecMech.Rev2PcetO2(k_rate=ks2, conc=[1.,1.])

    cyc1.pH = 1.
    cyc2.pH = 1.

    print(cyc1.rate(0.2))
    print(cyc2.rate(0.2))

    cyc_list = [cyc1, cyc2]

    mc = MultiCycle(cyc_list)

    mc.pH = 3

    print(mc.cyc_tuple[0].pH)
    print(mc.cyc_tuple[1].pH)

    #sim = SimTafel.SimTafel(mc) #not currently using any methods from this class

    rate_list = []
    for pH in pH_range:
        mc.pH = pH
        rate_list.append(np.log10(np.abs(mc.rate(V))))
        #cyc1.pH = pH
        #rate_list.append(np.log10(np.abs(cyc1.rate(V))))
        #rate_list.append(cyc1.rate(V))

    print(pH_range)
    print(rate_list)
    plt.plot(pH_range, rate_list)
    plt.show()
