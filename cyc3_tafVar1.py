import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import scipy.optimize as opt
import matplotlib as mpl

import cycleG

font = {'size':14}
mpl.rc('font',**font)

kT = 0.02585
h2o = -2.46
#G = np.array([0.,0.3,-1.8]) #free energy of intermediates
G = np.array([0.,-0.2,-1.8]) #free energy of intermediates
prod_G = np.array([0.,h2o,h2o]) #producing water changes dG for reaction step
dG = np.roll(G,-1) - G + prod_G #roll does a circular permuation of G
raw_barriers = np.array([0.1,0.01,0.1])

barriers = np.array([raw_barriers,raw_barriers])
for i,item in enumerate(dG):
    if np.sign(item) == 1: barriers[0,i] += item
    elif np.sign(item) == -1: barriers[1,i] -= item
A = np.ones((2,len(G)))*10.**8
k = A*np.exp(-1.*barriers/kT)

Vdep = np.array([1,0,1])
#pHdep = np.array([[0,0,0],[0,1,1]])
#pHdep = np.array([[0,0,0],[1,1,1]])
pHdep = np.array([[1,1,1],[0,0,0]])

cyc = cycleG.RcycG(k,Vdep,pHdep)

pH = 12
OH = 10.**(pH-14.)
V_range = np.linspace(-1.,1.5,200)
rlist = []
for V in V_range:
    rlist.append(cyc.rate(V,OH))
rlog = np.log10(rlist)


r_grad = np.gradient(rlog,V_range[1]-V_range[0])
fig, (ax0, ax1) = plt.subplots(nrows=2)
ax0.plot(V_range,rlog)
ax1.plot(V_range,1./r_grad)
ax1.set_ylim([-0.2,0])

plt.xlabel('V')

ax0.set_ylabel('$I$')
ax1.set_ylabel(r'log$_{10}$($I$)')

lw = 2
for ln in ax0.lines: ln.set_linewidth(lw)
for ln in ax1.lines: ln.set_linewidth(lw)
plt.show()









