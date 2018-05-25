"""
Half order pH dependence for mechanism fitting
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

import SimTafel
import ElecMech

dom = 100
pH_range = (3,14)
pH_dom = np.linspace(pH_range[0],pH_range[1],dom)

V_slope = 0.029 #dV(onset)/dpH at half order would be about this much. What mechanisms can explain this trend?
#V_slope = 0.059 #dV(onset)/dpH at half order would be about this much. What mechanisms can explain this trend?
onsetV_h = -1.*V_slope*pH_dom - 0.2



k_bounds = (0, 10**15)
start_pH = 3.

#k_param = (15,82.2)
#emech = ElecMech.BV(k_param, start_pH)

#k_param = (1.,2.,3.,4.)
#emech = ElecMech.RevPcet2(k_param,start_pH)

#k_param = (15,82.2,1.)
#emech = ElecMech.MultiBV(k_param,start_pH)

k_param = (1.,1.,1.,1.,1.,1.)
emech = ElecMech.RevPcet2ab(k_param,start_pH)




sim = SimTafel.SimTafel(emech)
popt, pcov = opt.curve_fit(sim.fitOnsetPH, pH_dom, onsetV_h, p0=k_param, bounds=k_bounds)    
print(popt)
emech.setConsts(popt)
fit_onset_Vs, grad_Vs = sim.onsetGradPH(pH_dom)


plt.plot(pH_dom, onsetV_h, label="29mV/pH Data")
plt.plot(pH_dom, fit_onset_Vs, label="Fit: Acid+Base BV")
plt.legend()
plt.xlabel("pH")
plt.ylabel(r"V$_{onset}$")
plt.show()
