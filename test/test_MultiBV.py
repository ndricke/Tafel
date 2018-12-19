import numpy as np
import matplotlib.pyplot as plt

import ElecMech
import SimTafel

dom = 16001
pH_dom = np.linspace(-3,14,dom)
V_dom = np.linspace(-2,2,dom)

mbv = ElecMech.MultiBV((10.**2,10.**2,1.), 9.) 

sim = SimTafel.SimTafel(mbv, ep=10**1)

#sim.plotTafel(V_dom)
sim.plotdVdPH(pH_dom)

"""Creates a figure that depicts the importance of the choice of j0"""
#fig, (ax0,ax1) = plt.subplots(nrows=2)
#
#for i in np.arange(-2,3):
#    sim = SimTafel.SimTafel(mbv, ep=10**i)
#    onset_Vs, grad_Vs = sim.onsetGradPH(pH_dom) 
#    ax0.plot(pH_dom, onset_Vs, label=r'j$_0$=10$^{%s}$' % str(i))
#    ax1.plot(pH_dom, grad_Vs)
#
#
#ax0.set_xlim([pH_dom[0],pH_dom[-1]])
#ax1.set_xlim([pH_dom[0],pH_dom[-1]])
#plt.xlabel('pH')
#ax0.set_ylabel(r'V')
#ax1.set_ylabel('$\partial$V/$\partial$pH')
#ax0.legend()


"""Creates scan of current vs voltage for varying pH to show transitions why V_onset has complex pH behavior"""
#fig, (ax0,ax1) = plt.subplots(nrows=2)
#
#for i in np.arange(-3,14):
#    mbv.pH = i
#    current = mbv.rate(V_dom)
#    #plt.plot(V_dom, np.log10(np.abs(current)), label="pH="+str(i))
#
#    logI = np.log10(np.abs(current))
#    dlogIdV = np.gradient(logI, V_dom[1]-V_dom[0]) #inverse of Tafel slope; calc this way bc V_dom indep var w/ even spacing
#    dVdlogI = 1./dlogIdV 
#
#    ax0.plot(V_dom, logI)
#    ax1.plot(V_dom, dVdlogI, label="pH="+str(i))
#    plt.xlabel('V')
#
##plt.ylim([-100,100])
##plt.xlim([-1,0])
##plt.legend()
##plt.ylabel("Current")
##plt.xlabel("Voltage")
#
#ax1.set_ylim([-0.15,0.15])
#ax0.set_xlim([V_dom[0], V_dom[-1]])
#ax1.set_xlim([V_dom[0], V_dom[-1]])
#ax0.set_ylabel('log(I)')
#ax1.set_ylabel('$\partial$V/$\partial$log(I)')
#ax1.legend(loc=4)

plt.show()
