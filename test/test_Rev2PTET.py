import numpy as np
import matplotlib.pyplot as plt

import ElecMech
import SimTafel

dom = 800
pH_dom = np.linspace(3,10,dom)
V_dom = np.linspace(-1,1,dom)

dVdPH_data = 2



#pcet2 = ElecMech.Rev2PTET((10.**8,10.**8,1.,1.,1.,10.**8), 1.) #Can create 60mV/dec and 120 mV/dec regimes in acidic conditions
pcet2 = ElecMech.Rev2PTET((10.**8,10.**-1,1.,1.,1.,10.**8), 14.) #Can we get only the 120 mV/dec regime to appear?

sim = SimTafel.SimTafel(pcet2)
sim.plotTafel(V_dom,[0,0.3] )
#sim.plotdVdPH(pH_dom)

print(sim.dVdlogI)



plt.show()






