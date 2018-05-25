import numpy as np
import matplotlib.pyplot as plt

import ElecMech
import SimTafel


dom = 400
pH_dom = np.linspace(0,18,dom)
V_dom = np.linspace(-3,3,dom)

cyc3 = ElecMech.Cyc3Trap((1.,1.,1.,1.,1.,1.), 3.) #Can create regime's with 120mV/dec transition to 40mV/dec

sim = SimTafel.SimTafel(cyc3, ep=10**-8)

#sim.plotTafel(V_dom)
sim.plotdVdPH(pH_dom)


plt.show()

