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
k_list = [1,0.01,0.01,1] 
dE_range = [-4,4]
sig = 0.01

dis = IH.DisorderGCC(k_list, sig)
dis.pH = 6
sig_list = [0.01, 0.05, 0.1, 0.2, 0.4]

sim = SimTafel.SimTafel(dis)

## pH dependence with intrinsic disorder
"""
V_list = []
onset_J = 0.1
n = 20
logO2_list = np.linspace(-12, -1, n)

fig, (ax0,ax1) = plt.subplots(nrows=2)

for sig in sig_list:
    dis.sig = sig
    V_arr = np.zeros(n)
    for i, logO2 in enumerate(logO2_list):
        O2 = 10**(logO2)
        dis.O2 = O2
        V_arr[i] = sim.findOnsetV(onset_J, onsetV_guess=-0.03)['x']

    ax0.plot(logO2_list, V_arr, label=str(sig))
    ax0.scatter(logO2_list, V_arr)

    dVdO2 = np.gradient(V_arr, logO2_list[1]-logO2_list[0]) #inverse of Tafel slope; calc this way bc V_dom indep var w/ even spacing
    ax1.plot(logO2_list, dVdO2)
    ax1.scatter(logO2_list, dVdO2)

plt.xlabel("log(O2)")
ax0.set_ylabel("Onset Potential (V)")
ax1.set_ylabel("dVdO2")
#plt.xlim([-0.01, 0.21])
fig.set_size_inches(11.,11.,forward=True)
ax0.legend(loc=2)
#plt.savefig("DisorderGCC_1t1000_pHvsOnsetV.png", transparent=True, bbox_inches='tight', pad_inches=0.02)
plt.show()
"""

## LogO2 vs Rate

V = -0.6
n = 20
logO2_list = np.linspace(-12, 2, n)
fig, (ax0,ax1) = plt.subplots(nrows=2)

for sig in sig_list:
    dis.sig = sig
    rate_arr = np.zeros(n)
    for i, logO2 in enumerate(logO2_list):
        O2 = 10**(logO2)
        dis.O2 = O2
        rate_arr[i] = dis.rate(V)

    ax0.plot(logO2_list, np.log10(rate_arr))
    ax0.scatter(logO2_list, np.log10(rate_arr))

    dRdO2 = np.gradient(np.log10(rate_arr), logO2_list[1]-logO2_list[0])
    ax1.plot(logO2_list, dRdO2, label=str(sig))
    ax1.scatter(logO2_list, dRdO2)

plt.xlabel(r"log(O$_2$)")
ax0.set_ylabel("log(Rate)")
ax1.set_ylabel("d(log(Rate))/d(log(O2))")
ax1.legend(loc=1)
plt.show()





