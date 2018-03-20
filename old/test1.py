import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import scipy.optimize as opt

f = 38.684
V = np.linspace(-1,1,400)
r = np.exp(-0.5*f*V)
rlog = np.log10(r)
rgrad = np.gradient(rlog,V[1]-V[0])

#plt.plot(V,rlog)
plt.plot(V,1./rgrad)

#plt.plot(rlog,V)
#plt.plot(V,V*-0.059)


x = np.linspace(-1,1,1000)
y = -x/0.118
grady = np.gradient(y,x[1]-x[0])
#plt.plot(x,y)
plt.plot(x,1./grady)

plt.show()
