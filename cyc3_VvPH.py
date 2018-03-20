import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy as sc
import scipy.optimize as opt

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

V = 0.
pH_range = [-2.,8]
h = np.linspace(pH_range[0],pH_range[1],200)
H = 10.**(-1.*h)
OH = 10.**(h-14.)

V_list = [0] #start with a 0 so the code can cleanly use this as an initial guess
for i in H:
    cyc.H = i
    V_sol = opt.root(cyc.varV,V_list[-1],method='hybr')
    V_list.append(V_sol.x[0])
V_list.pop(0)
gradV = np.gradient(V_list,h[1]-h[0])
print gradV

#plt.plot(h,V_list)
fig, (ax0, ax1) = plt.subplots(nrows=2)
ax0.plot(h,V_list)
ax1.plot(h,gradV)
ax0.set_ylabel(r'V')
ax1.set_ylabel('$\partial$V/$\partial$pH')
lw = 2
for ln in ax0.lines: ln.set_linewidth(lw)
for ln in ax1.lines: ln.set_linewidth(lw)

plt.xlabel('pH')
ax0.set_ylabel(r'V')
ax1.set_ylabel('$\partial$V/$\partial$pH')
plt.show()




























