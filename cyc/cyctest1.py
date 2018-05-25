import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy as sc
import scipy.optimize as opt

import cycleG

h2o = -2.46

#G = np.array([0.,-0.1,-0.1]) #free energy of intermediates
#prod_G = np.array([0.,h2o,h2o]) #producing water changes dG for reaction step
#raw_barriers = np.array([0.5,0.01,0.01])
#A = 10.**8
#Vdep = np.array([1,1,1])
#pHdep = np.array([[1,1,1],[0,0,0]])

G = np.array([0.,0.]) #free energy of intermediates
prod_G = np.array([0.,h2o]) #producing water changes dG for reaction step
raw_barriers = np.array([0.5,0.01])
A = 10.**8
Vdep = np.array([1,1])
pHdep = np.array([[1,1],[0,0]])

k = cycleG.genK(G,prod_G,raw_barriers,A)
cyc = cycleG.RcycG(k,Vdep,pHdep)














