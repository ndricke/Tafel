import scipy as sc
import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt

import ElecMech
import SimTafel

#fun = lambda x: (x[0] - 1)**2 + (x[1] - 2.5)**2
#cons = ({'type': 'ineq', 'fun': lambda x: x[0] - 2 * x[1] + 2})
#bnds = ((0, None), (0,None))
#res = opt.minimize(fun, (2,0), method='SLSQP', bounds=bnds, constraints=cons)
#print(res.x)

#def parab(x, a, b, c):
#    return a*x**2 + b*x + c
#
#def cub(x, a, b, c, d):
#    return a*x**3 + b*x**2 + c*x + d
#
#class fitter:
#
#    def __init__(self, domain, data):
#        self.domain = domain
#        self.data = data
#
#    def fit_data(self, p):
#        a, b, c = p
#        return np.linalg.norm(self.data - parab(self.domain, a, b, c))
#
#
#dom = np.linspace(1,7,500)
#cubix = cub(dom, 0.5, -2.8, 1.11, 3.23)
#
#fit = fitter(dom, cubix)
#res = opt.minimize(fit.fit_data, x0=(-2.8, 1.11, 3.23)) #, method='SLSQP')
#print(res.x)
#a_opt, b_opt, c_opt = res.x
#
#plt.plot(dom, parab(dom, a_opt, b_opt, c_opt))
#plt.plot(dom, cubix)
#plt.show()

class DataFit:
    def __init__(s, Vdom, Jdata, emech):
        s.emech = emech
        s.Jdata = Jdata
        s.Vdom = Vdom

    def evalfit(s, kconsts):
        s.emech.setConsts(kconsts)
        return np.linalg.norm(s.emech.rate(s.Vdom) - s.Jdata)

dom = 800

data = np.loadtxt('echem_data.csv',delimiter=',') #data[:,0] = log10(J), data[:,1] = V
data[:,1] *= -1.

J10_data = 10.**data[:,0]
V_dom = np.linspace(data[0,1], data[-1,1], dom)
V_dom_exp = np.linspace(-2,2, dom)

#r1 = 10.**-3
#r2 = 10.**-7*2.7
#r_tup = (r1,r2)
#bnds = ((0.,None), (0.,None))
#mech = ElecMech.ClassicBV((r1,r2), 0.)

r_tup = (10., 0.1, -0.5, 0.8, 1.2)
bnds = ((0, None),(-2,2),(-2,2),(0,20),(0,20))
mech = ElecMech.RevPcet2MHC(r_tup, 1.) 

#r_tup = (0.1, -0.4, 0.8)
#bnds = ((0, None),(None,None),(0,None))
#mech = ElecMech.MHC(r_tup, 1.) 

fit = DataFit(Vdom=data[:,1], Jdata=J10_data, emech=mech)
res = opt.minimize(fit.evalfit, x0=r_tup, bounds=bnds, method='SLSQP')
#res = opt.minimize(fit.evalfit, x0=r_tup, bounds=bnds, method='TNC')

print(res.x)
print(res.success)
mech.setConsts(res.x)
print(fit.evalfit(res.x))


#plt.plot(V_dom, mech.rate(V_dom))
plt.plot(V_dom_exp, mech.rate(V_dom_exp))

plt.scatter(data[:,1], J10_data)
plt.xlabel('V')
plt.ylabel('J')

plt.show()


























