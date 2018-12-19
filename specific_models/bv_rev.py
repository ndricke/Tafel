from BV import BV

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy as sc
import scipy.optimize as opt

matplotlib.rc('xtick', labelsize=20) 
matplotlib.rc('ytick', labelsize=20) 


#font = {'family' : 'normal',
#        'weight' : 'bold',
#        'size'   : 22}
#matplotlib.rc('font', **font)

Vs = np.linspace(-2,-0.3,400)

cyc = BV(a=0.01, b=1.)

"""Create Tafel plots over a range of pH"""
#pH_marks = np.arange(3,14,2)
#Is_pH = []
#for pH in pH_marks:
#    cyc.pH = pH
#    Is_pH.append(cyc.rate(Vs))
##    plt.plot(Vs, Is_pH[-1])
#    plt.plot(np.log10(Is_pH[-1]), Vs, label="pH = "+str(pH), linewidth=2.0)
#
##plt.ylim([-5,5])
#plt.xlabel(r'Log$_(10)$(I)', fontsize=20)
#plt.ylabel('V', fontsize=20)
#plt.legend()
#plt.show()


"""Create a plot showing how the pH dependence changes"""
V = 0.
pH = np.linspace(3.,14.,200)
H = 10.**(-pH)

V_onset = cyc.fitOnset(pH, 1, 1)
gradV = np.gradient(V_onset,pH[1]-pH[0])

fig, (ax0,ax1) = plt.subplots(nrows=2)
ax0.plot(pH,V_onset)
ax1.plot(pH,gradV)

print(gradV)

ax1.set_ylim([-0.2,0.])

ax0.set_xlim([3,14])
ax1.set_xlim([3,14])
plt.xlabel('pH', fontsize=20)
ax0.set_ylabel(r'V', fontsize=20)
ax1.set_ylabel('$\partial$V/$\partial$pH', fontsize=20)
lw = 2
for ln in ax0.lines: ln.set_linewidth(lw)
for ln in ax1.lines: ln.set_linewidth(lw)
plt.show()

