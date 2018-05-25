import numpy as np
import matplotlib.pyplot as plt

import ElecMech
import SimTafel

dom = 16001
V_dom = np.linspace(-0.2,0.2,dom)

bv = ElecMech.ClassicBV((10.**12,10.**12)) 
sim = SimTafel.SimTafel(bv, ep=10**1)

fig, (ax0,ax1) = plt.subplots(nrows=2)
for i in range(12):
    bv.setConsts((10.**i,10.**i))
    current = bv.rate(V_dom)

    logI = np.log10(np.abs(current))
    dlogIdV = np.gradient(logI, V_dom[1]-V_dom[0]) #inverse of Tafel slope; calc this way bc V_dom indep var w/ even spacing
    dVdlogI = 1./dlogIdV 

    ax0.plot(V_dom, logI)
    ax1.plot(V_dom, dVdlogI, label="k=10^"+str(i))
    plt.xlabel('V')

ax0.set_ylabel('log(I)')
ax1.set_ylabel('$\partial$V/$\partial$log(I)')
ax1.legend(loc=4)

plt.show()


