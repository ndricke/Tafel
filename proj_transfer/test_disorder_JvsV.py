"""
This script shows that disorder can affect the rate of a reaction without affecting the tafel slope
The specific reaction is HER, where this model calculations the current-voltage relationships for varying disorder
Since both reaction steps have the same current-voltage relationship, mixing them does not affect the slopes
"""

import matplotlib.pyplot as plt
import numpy as np
from DisorderMech import *


V_dom = np.linspace(-0.2,0.2,500) # domain for applied voltage
dE_range = [-4,4]
sig_list = [0.01, 0.05, 0.1, 0.2, 0.4] # list of values for varying disorders

## k1, kn1, k2, kn2
k_list = [1,1000,1000,1]

dis = DisorderHER(k_list, sig_list[0])
rate_list = []
disorder_list = []
dis.H = 1.

sim = SimTafel.SimTafel(dis)

## Plot Tafel slope as a function of intrinsic disorder

fig, (ax0,ax1) = plt.subplots(nrows=2)
for s in sig_list:
    dis.sig = s
    logI_list = []; dV_list = [];
    for V in V_dom:
        I = dis.rate(V)
        logI = np.log10(np.abs(I))
        logI_list.append(logI)

    dlogIdV = np.gradient(logI_list, V_dom[1]-V_dom[0]) #inverse of Tafel slope; calc this way bc V_dom indep var w/ even spacing
    dVdlogI = 1./np.array(dlogIdV)

    ax0.plot(V_dom, logI_list)
    ax1.plot(V_dom, dVdlogI, label=str(s))

plt.xlabel('V')
#ax1.set_ylim(input_ylim)
ax0.set_xlim([V_dom[0], V_dom[-1]])
ax1.set_xlim([V_dom[0], V_dom[-1]])
ax0.set_ylabel('log(I)')
ax1.set_ylabel('$\partial$V/$\partial$log(I)')
ax1.legend()
fig.set_size_inches(11.,11., forward=True)

plt.show()
#plt.savefig("IntrinsicHet_1t1_VvsJ.png", transparent=True, bbox_inches='tight', pad_inches=0.02)
