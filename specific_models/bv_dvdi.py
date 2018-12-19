import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import scipy.optimize as opt

from cycleG import BV

V = np.linspace(-0.8,-0.11,400)
H = 1.

cyc = BV(10**-6)
I_list = []
for i in V:
    I_list.append(cyc.rate(V=i))

I = np.array(I_list)
logI = np.log10(I)

#dIdV = np.gradient(I,V[1]-V[0])
#dVdI = 1./dIdV

dlogIdV = np.gradient(logI,V[1]-V[0])
dVdlogI = 1./dlogIdV


#plt.plot(V,I)

fig, (ax0,ax1) = plt.subplots(nrows=2)
ax0.plot(V,I)
#ax1.plot(V,dVdI)
#ax1.plot(V,dIdV)
ax1.plot(V,dVdlogI)
#ax1.set_ylim([-0.2,0.])
#ax0.set_xlim([3,14])
#ax1.set_xlim([3,14])
plt.xlabel('V')
ax0.set_ylabel(r'I')
ax1.set_ylabel('$\partial$V/$\partial$log10(I)')
lw = 2
for ln in ax0.lines: ln.set_linewidth(lw)
for ln in ax1.lines: ln.set_linewidth(lw)
plt.show()

