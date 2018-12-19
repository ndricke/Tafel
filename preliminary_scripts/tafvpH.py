import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import scipy.optimize as opt


def tvp(A,h,E):
    a = 0.5
    f = 38.949
    return (np.log(A*10.**-h+1)-np.log(E))/(a*f)

def tvp2(A,h):
    return np.log(A*10.**-h+1)

class TvpRev(object):
    a = 0.5
    f = 38.949
    ep = 10**-1
    def __init__(s,h,A,B,th1,th2):
        s.h = h
        s.A = A
        s.B = B
        s.th1 = th1
        s.th2 = th2

    def overp(s, eta):
        fwd = (s.A*10.**-s.h + s.B)*np.exp(-s.a*eta*s.f)*s.th1
        bck = (s.A + s.B*10.**s.h)*np.exp((1.-s.a)*eta*s.f)*s.th2
        return fwd - bck - s.ep
  

#A = 100
#h = np.linspace(0,4,50)
#eta = tvp(A,h,10**-4)
#plt.plot(h,eta)
#plt.show()

pH_range = [-9,14]
h = np.linspace(pH_range[0],pH_range[1],200)
#A = 1
#B = 1
A = 0.001
B = 0.00001
th1 = 0.1
th2 = 0.9
sol_list = [0] #start with a 0 so the code can cleanly use this as an initial guess
for i in h:
    pr = TvpRev(i,A,B,th1,th2)
#    sol = opt.root(pr.overp,sol_list[-1],method='lm')
    sol = opt.root(pr.overp,sol_list[-1],method='hybr')
    sol_list.append(sol.x[0])
sol_list.pop(0)

print sol_list
plt.plot(h,sol_list)
plt.xlabel("Relative pH")
plt.ylabel("Overpotential (V)") #when current is observable
plt.xlim(pH_range)
plt.show()

#n = np.linspace(-0.2,0.2,200)
#pr = TvpRev(0.2,A,B,th1,th2)
#n_scan = pr.overp(n)
#plt.plot(n,n_scan)
#plt.show()
















