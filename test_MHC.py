import numpy as np
import matplotlib.pyplot as plt

import ElecMech
import SimTafel


f = 38.9



dom = 4000
pH_dom = np.linspace(-8,18,dom)
V_dom = np.linspace(-4,4,dom)

mhc = ElecMech.MHC((1.,2.3,2.3), 0.)
#mhc = ElecMech.RevPcet2MHC((10/f, -20/f, 55./f, 35./f), 0.) #pretty interesting, but the regimes we care about are small
#mhc = ElecMech.RevPcet2MHC((1.0, 0.1, -1.0, 0.5, 1.), 0.) 

sim = SimTafel.SimTafel(mhc)


I = mhc.rate(V_dom)
logI = np.log10(np.abs(I))

dlogIdV = np.gradient(logI, V_dom[1]-V_dom[0]) #inverse of Tafel slope; calc this way bc V_dom indep var w/ even spacing
dVdlogI = 1./dlogIdV #A flat 120 mV/dec tafel slope for negative V's

fig, (ax0,ax1) = plt.subplots(nrows=2)
ax0.plot(V_dom, logI)
ax1.plot(V_dom, dVdlogI)
plt.xlabel('V')
ax1.set_ylim([-.2,.2])

#plt.plot(V_dom, I)

plt.show()



