import numpy as np
import matplotlib.pyplot as plt

import ElecMech
import SimTafel






dom = 4000
pH_dom = np.linspace(-8,18,dom)
V_dom = np.linspace(0.,0.1,dom)

pcet2 = ElecMech.RevPcet2((1.,1.,1.,1.), 0.) #Can create regime's with 120mV/dec transition to 40mV/dec

sim = SimTafel.SimTafel(pcet2)

#onset_Vs, grad_Vs = sim.OnsetGradPH(pH_dom) #Can create regimes with 59 mV/pH, 80 mV/pH, and 120 mV/pH
#
#fig, (ax0,ax1) = plt.subplots(nrows=2)
#ax0.plot(pH_dom, onset_Vs)
#ax1.plot(pH_dom, grad_Vs)


#I = pcet2.rate(V_dom)
#logI = np.log10(np.abs(I))
#
#dlogIdV = np.gradient(logI, V_dom[1]-V_dom[0]) #inverse of Tafel slope; calc this way bc V_dom indep var w/ even spacing
#dVdlogI = 1./dlogIdV #A flat 120 mV/dec tafel slope for negative V's
#
#fig, (ax0,ax1) = plt.subplots(nrows=2)
#ax0.plot(V_dom, logI)
#ax1.plot(V_dom, dVdlogI)
#plt.xlabel('V')
#ax1.set_ylim([-0.2,0.2])

sim.plotTafel(V_dom, input_ylim=[0.,0.15])
#sim.plotdVdPH(pH_dom)


plt.show()




