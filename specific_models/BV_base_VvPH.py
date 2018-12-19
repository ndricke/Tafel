import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import scipy.optimize as opt

class BV(object):
    f = 38.684
    def __init__(s,ep): 
        s.ep = ep #what is the smallest point where they can observe the current?

    def rate(s,V,H=1.,f=38.684):
        return 0.01*H*np.exp(-0.5*f*V) - np.exp(0.5*f*V)
#        return H*np.exp(-0.5*f*V)

    def rate2(s,V,H=1.,f=38.684):
        OH = (10.**-14)/H
        return 0.00001*np.exp(-0.5*f*V) - OH*np.exp(0.5*f*V)

    def varV(s,V):
        return s.rate2(V,s.H) - s.ep

#V = np.linspace(-1,1,400)
#H = 1.

V = 0.
pH = np.linspace(0.,14.,200)
H = 10.**(-pH)

cyc = BV(10**-6)
V_list = [0] #start with a 0 so the code can cleanly use this as an initial guess
for i in H:
    cyc.H = i
    V_sol = opt.root(cyc.varV,V_list[-1],method='hybr')
    V_list.append(V_sol.x[0])
V_list.pop(0)

gradV = np.gradient(V_list,pH[1]-pH[0])

fig, (ax0,ax1) = plt.subplots(nrows=2)
ax0.plot(pH,V_list)
ax1.plot(pH,gradV)

ax1.set_ylim([-0.2,0.])

ax0.set_xlim([0,14])
ax1.set_xlim([0,14])
plt.xlabel('pH')
ax0.set_ylabel(r'V')
ax1.set_ylabel('$\partial$V/$\partial$pH')
lw = 2
for ln in ax0.lines: ln.set_linewidth(lw)
for ln in ax1.lines: ln.set_linewidth(lw)
plt.show()

