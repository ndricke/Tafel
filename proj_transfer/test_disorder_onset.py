"""
This plot shows the onset potential as a function of disorder
If the rate constants are all at optimal values, namely so that the reaction is reversible, disorder will not increase the rate
If the rate constants are suboptimal, an even spreading of the rate constants will increase the rate (up to a point)
"""

import matplotlib.pyplot as plt
import numpy as np
from DisorderMech import *

## pH dependence with intrinsic disorder

V_dom = np.linspace(-0.2,0.2,500) # domain for applied voltage
dE_range = [-4,4]
sig = 0.001

## k1, kn1, k2, kn2
k_list = [1,1000,1000,1]

dis = DisorderHER(k_list, sig)
rate_list = []
disorder_list = []
dis.H = 1.

sim = SimTafel.SimTafel(dis)



## Plot onset potential as a function of intrinsic disorder
V_list = []
onset_J = 0.1
n = 60
sig_list = np.linspace(0.001, 0.2, n)

fig, ax = plt.subplots()

## Test to make sure that disorder converges to ordered at sig = 0.
#ordered_mech = ElecMech.RevPcet2()
#ordered_mech.setConsts(k_list)
#ordered_mech.setConc(0)
#ordered_sim = SimTafel.SimTafel(ordered_mech)
#V_ordered = ordered_sim.findOnsetV(onset_J)['x']
#plt.scatter(0., V_ordered)

k_ratios = [1, 5, 10, 25, 50, 100]
for ratio in k_ratios:
    ks = [1./ratio,ratio,ratio,1./ratio]
    dis.setConsts(ks)
    print(dis.k1, dis.kn1, dis.k2, dis.kn2)
    V_arr = np.zeros(n)
    for i, sig in enumerate(sig_list):
        dis.sig = sig
        V_arr[i] = sim.findOnsetV(onset_J, onsetV_guess=-0.03)['x']

    plt.scatter(sig_list, V_arr, )
    plt.plot(sig_list, V_arr, label=ratio)

plt.xlabel("Disorder")
plt.ylabel("Onset Potential (V)")
plt.xlim([-0.01, 0.21])
fig.set_size_inches(11.,11.,forward=True)
plt.legend()
#plt.savefig("IntrinsicHet_DisorderVsOnsetV.png", transparent=True, bbox_inches='tight', pad_inches=0.02)
plt.show()
