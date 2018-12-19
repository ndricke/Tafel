import numpy as np
import matplotlib.pyplot as plt

import ElecMech
import SimTafel






dom = 400
pH_dom = np.linspace(-15,30,dom)
V_dom = np.linspace(-1,1,dom)

mech = ElecMech.RevPcet2ab(conc=3.) #Can create regime's with 120mV/dec transition to 40mV/dec
mech.genConsts((-0.5, 0.2, 0.3, 0.1, 0.4, 0.2))


sim = SimTafel.SimTafel(mech, ep=10**-6)

#sim.plotTafel(V_dom)
sim.plotdVdPH(pH_dom)


plt.show()




