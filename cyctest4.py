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

#gG = np.array([0.,0.3,-1.8]) #free energy of intermediates
#G = np.array([0.,-0.2,-1.8]) #free energy of intermediates
#prod_G = np.array([0.,h2o,h2o]) #producing water changes dG for reaction step
#dG = np.roll(G,-1) - G + prod_G #roll does a circular permuation of G
#raw_barriers = np.array([0.1,0.01,0.1])

G = np.array([0.,-0.2,-1.8]) #free energy of intermediates
prod_G = np.array([0.,h2o,h2o]) #producing water changes dG for reaction step
dG = np.roll(G,-1) - G + prod_G #roll does a circular permuation of G
raw_barriers = np.array([0.1,0.01,0.1])
A = 10.**8

#barriers = np.array([raw_barriers,raw_barriers])
#for i,item in enumerate(dG):
#    if np.sign(item) == 1: barriers[0,i] += item
#    elif np.sign(item) == -1: barriers[1,i] -= item
#A = np.ones((2,len(G)))*10.**8
#k = A*np.exp(-1.*barriers/kT)

k = cycleG.genK(G,prod_G,raw_barriers,A)
print k

Vdep = np.array([1,1,1])
#pHdep = np.array([[0,0,0],[0,1,1]])
#pHdep = np.array([[0,0,0],[1,1,1]])
pHdep = np.array([[1,1,1],[0,0,0]])

cyc = cycleG.RcycG(k,Vdep,pHdep)

pH_range = [-2.,14.]
h = np.linspace(pH_range[0],pH_range[1],200)
H = 10.**(-1.*h)

V_list = [0] #start with a 0 so the code can cleanly use this as an initial guess
for i in H:
    cyc.H = i
    V_sol = opt.root(cyc.varV,V_list[-1],method='hybr')
    V_list.append(V_sol.x[0])
V_list.pop(0)

#print sol_list

fig, (ax0,ax1) = plt.subplots(nrows=2)
ax0.plot(h,V_list)
grad_V = np.gradient(V_list,h[1]-h[0])
#print grad_V
ax1.plot(h,grad_V)
plt.xlabel("pH")
plt.xlim(pH_range)

ax0.set_ylabel("Potential (V)") #when current is observable
ax1.set_ylabel("$\partial$V/$\partial$pH") #when current is observable

lw = 2
for ln in ax0.lines: ln.set_linewidth(lw)
for ln in ax1.lines: ln.set_linewidth(lw)
plt.show()



























