import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

import ElecMech
import SimTafel


class DataFit:
    def __init__(s, data, emech):
        s.emech = emech
        s.data = data

    def evalfit(s, kconsts):
        s.emech.setConsts(kconsts)
        return np.linalg.norm(s.emech.lograte(data[:,1]) - data[:,0])



dom = 800

data = np.loadtxt('echem_data.csv',delimiter=',') #data[:,0] = log10(J), data[:,1] = V
V_dom = np.linspace(data[0,1], data[-1,1], dom)
#V_dom = np.linspace(-2,2, dom)

r1 = 10.**-3
r2 = 10.**-7*2.7
con = 10.**4
r1, r2 = r1*con, r2*con
#mech = ElecMech.ClassicBV((r1,r2), 0.)


#r1, r2, r3, r4 = 
#mech = ElecMech.RevPcet2((r1, r1, r2, r2), 1.) 

k_list = (100., 0.2, -0.5, 2.3, 1.2)
k_bounds = ((0, 10**10),(-2,2),(-2,2),(0,20),(0,20))
mech = ElecMech.RevPcet2MHC(k_list, 1.) 

#print(V_dom)
#print(mech.rate(V_dom))
#print(mech.lograte(V_dom))

fitter = DataFit(data, mech)
#r = [r1, r2]
#print(fitter.evalfit(r))



#sim = SimTafel.SimTafel(mech)
#popt, pcov = opt.curve_fit(sim.fitLograte, data[1], data[0], p0=(10.**1, 10**-3*2.7), bounds=k_bounds)    
#print(popt)
#mech.setConsts(popt)



#res = opt.minimize( fitter.evalfit, x0=[10.**1,10.**-3], bounds=k_bounds)    
#res = opt.minimize( fitter.evalfit, x0=[10.**1,10.**-3])    
#res = opt.minimize( fitter.evalfit, x0=[r1,r1,r2,r2])
#res = opt.minimize( fitter.evalfit, x0=[0.2, -0.5, 2.3, 1.2])
res = opt.minimize( fitter.evalfit, x0=k_list, bounds=k_bounds)

#res = opt.minimize( fitter.evalfit, x0=[10.**1, 2.7*10.**-3])

#mech.setConsts(res.x)
#print(res.x)
#print(fitter.evalfit(res.x))





#def rosen(x):
#    return sum(100.0*(x[1:] - x[:-1])**2.0 + (1 - x[:-1]**2.0)**2)
#x0 = [1.3, 0.7, 0.8, 1.9, 1.2]
#res = opt.minimize(rosen, x0)
print(res.x)

plt.scatter(data[:,0], data[:,1])
plt.plot(mech.lograte(V_dom), V_dom)
plt.ylabel('V')
plt.xlabel('log(J)')
plt.show()






