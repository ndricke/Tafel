import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import scipy.optimize as opt

import cycleG

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
rlist = []
#V_range = np.linspace(-1.,1.,200)
V = 0.

#for V in V_range:
#    rlist.append(cyc.rate(V,OH))
#rlog = np.log10(rlist)

pH_range = [2.,14]
h = np.linspace(pH_range[0],pH_range[1],200)
H = 10.**(-1.*h)
OH = 10.**(h-14.)

#Plot pH dependence
for H_ind in H:
    rlist.append(cyc.rate(V,H_ind))
plt.plot(h,np.log10(rlist))

rlist = []; V = -0.2
for H_ind in H:
    rlist.append(cyc.rate(V,H_ind))
plt.plot(h,np.log10(rlist))

rlist = []; V = -0.4
for H_ind in H:
    rlist.append(cyc.rate(V,H_ind))
plt.plot(h,np.log10(rlist))

rlist = []; V = -0.6
for H_ind in H:
    rlist.append(cyc.rate(V,H_ind))
plt.plot(h,np.log10(rlist))

#sol_list = [0] #start with a 0 so the code can cleanly use this as an initial guess
#for i in OH:
#    cyc.OH = i
#    sol = opt.root(cyc.varV,sol_list[-1],method='hybr')
#    sol_list.append(sol.x[0])
#sol_list.pop(0)

#sol_list = [0] #start with a 0 so the code can cleanly use this as an initial guess
#for i in H:
#    cyc.H = i
#    sol = opt.root(cyc.varV,sol_list[-1],method='hybr')
#    sol_list.append(sol.x[0])
#sol_list.pop(0)
#
#plt.plot(h,sol_list)
##plt.plot(h,np.gradient(sol_list))
#plt.xlabel("pH")
#plt.ylabel("Potential (V)") #when current is observable
#plt.xlim(pH_range)

#x = np.linspace(4,14,200)
#plt.plot(x,-x+5)
#plt.plot(x,-2.*x+5)

plt.show()




























