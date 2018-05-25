
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

import SimTafel
import ElecMech




dom = 100
pH_range = (3,14)
pH_dom = np.linspace(pH_range[0],pH_range[1],dom)

V_slope = 0.059 
onsetV_h = -1.*V_slope*pH_dom - 0.2


mbv = ElecMech.MultiBV((15,82.2,1.),3)
sim = SimTafel.SimTafel(mbv)
popt, pcov = opt.curve_fit(sim.fitOnsetPH, pH_dom, onsetV_h, p0=(1.,1.,1.), bounds=(0.,10**12))    
print(popt)
mbv2 = ElecMech.MultiBV((popt[0],popt[1],popt[2]),3)
simBV2 = SimTafel.SimTafel(mbv2)
fit_onset_Vs, grad_Vs = simBV2.onsetGradPH(pH_dom)

















