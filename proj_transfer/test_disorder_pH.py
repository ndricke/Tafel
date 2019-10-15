"""
pH dependence of onset potential with intrinsic disorder
The specific reaction is HER
The plot this generates shows that disorder can change the point at which the pH dependence changes,
but it does not significantly influence this transition point for HER with the given parameters
"""

import matplotlib.pyplot as plt
import numpy as np
from DisorderMech import *

## pH dependence with intrinsic disorder

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

V_list = []
onset_J = 0.1
n = 40
pH_list = np.linspace(-2, 8, n)

fig, (ax0,ax1) = plt.subplots(nrows=2)

for sig in sig_list:
    dis.sig = sig
    V_arr = np.zeros(n)
    for i, pH in enumerate(pH_list):
        dis.pH = pH
        V_arr[i] = sim.findOnsetV(onset_J, onsetV_guess=-0.03)['x']

    ax0.plot(pH_list, V_arr)

    dVdpH = np.gradient(V_arr, pH_list[1]-pH_list[0]) #inverse of Tafel slope; calc this way bc V_dom indep var w/ even spacing
    ax1.plot(pH_list, dVdpH, label=str(sig))

plt.xlabel("pH")
ax0.set_ylabel("Onset Potential (V)")
ax1.set_ylabel("dV/dpH")
#plt.xlim([-0.01, 0.21])
fig.set_size_inches(11.,11.,forward=True)
plt.legend()
#plt.savefig("IntrinsicHet_1t1000_pHvsOnsetV.png", transparent=True, bbox_inches='tight', pad_inches=0.02)
plt.show()
