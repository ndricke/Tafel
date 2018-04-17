"""
Half order pH dependence for mechanism fitting
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

import simTafel

dom = 400
pH_range = (3,14)
pH_dom = np.linspace(pH_range[0],pH_range[1],dom)

V_slope = 0.029 #dV(onset)/dpH at half order would be about this much. What mechanisms can explain this trend?


onsetV_h = -1.*V_slope*pH_dom - 0.2


bv = simTafel.BV(a=15,b=82.2)
popt, pcov = opt.curve_fit(bv.fitOnset, pH_dom, onsetV_h, p0=(1.,1.))    
print(popt)
bv2 = simTafel.BV(a=popt[0], b=popt[1])
simBV2 = simTafel.SimTafel(bv2)
fit_onset_Vs, grad_Vs = simBV2.OnsetScanPH(pH_dom)

#rbv = simTafel.MultiBV()
#popt, pcov = opt.curve_fit(rbv.fitOnset, pH_dom, onsetV_h, p0=(1.,1.,1.))    
#print(popt)
#print(rbv.A, rbv.B, rbv.C)

##rbv2 = simTafel.MultiBV.fitOnset(pH_dom, A=popt[0], B=popt[1], C=popt[2])
##rbv2 = simTafel.MultiBV.fitOnset(pH_dom, (popt[0], popt[1], popt[2]))
#sim_rbv = simTafel.SimTafel(rbv)
#fit_onset_Vs, grad_Vs = sim_rbv.OnsetScanPH(pH_dom)



plt.plot(pH_dom, onsetV_h)
plt.plot(pH_dom, fit_onset_Vs)
plt.show()
