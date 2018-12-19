import sys

import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.integrate as integrate

import intrinsic_het as IH
import SimTafel

V_dom = np.linspace(-0.2,0.2,500)
## k1, kn1, k2, kn2
k_list = [1,1000,1000,1] 
dE_range = [-4,4]
sig = 0.01

dis = IH.DisorderPEET(k_list, sig)
sig_list = [0.01, 0.05, 0.1, 0.2, 0.4,0.5,0.6]

sim = SimTafel.SimTafel(dis)

## pH dependence with intrinsic disorder
V_list = []
onset_J = 0.1
n = 80
pH_list = np.linspace(-2, 12, n)

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
#plt.savefig("DisorderPEET_1t1000_pHvsOnsetV.png", transparent=True, bbox_inches='tight', pad_inches=0.02)
plt.show()





