import numpy as np
import matplotlib.pyplot as plt


#the two possible paths have the same free energy of reaction, but may have different barriers (and prefactors?)
#we can absorb the dG between the states into the overpotential and pH dependence
#explicitly writing the pH dependence should give us Nernstian behavior, and we are looking for deviations from it
#in particular, we would expect the deviation to occur while one rate takes over for the other
#Activation energy can be absorbed into the prefactor
class ARevPH:
    def __init__(s,A_H,A_OH,a=0.5):
        s.A_H = A_H
        s.A_OH = A_OH
        s.a = a
        s.f = 38.949

#pH just appears as a shift in overpotential; for eta=0 and pH = 0 the free energy difference bw states is 0
    def rate(s,eta,pH,cov):
        ka = s.A_H*np.exp(-(0.059*pH+s.a*eta)*s.f)
        kna = s.A_H*np.exp((1.0-s.a)*eta*s.f)
        kb = s.A_OH*np.exp(-s.a*eta*s.f)
        knb = s.A_OH*np.exp((0.059*pH+(1-s.a)*eta)*s.f)
        return (ka + kb)*cov - (kna + knb)*(1.-cov),ka,kna,kb,knb

class RevPH:
    def __init__(s,ka,kna,mech_ratio,a=0.5):
        s.koa = ka
        s.kona = kna
        s.kob = ka*mech_ratio
        s.konb = kna*mech_ratio
        s.a = a
        s.f = 38.949

    def rate(s,eta,pH,cov):
        ka = s.koa*np.exp(-(0.059*pH+s.a*eta)*s.f)
        kna = s.kona*np.exp((1.0-s.a)*eta*s.f)
        kb = s.kob*np.exp(-s.a*eta*s.f)
        knb = s.konb*np.exp((0.059*pH+(1-s.a)*eta)*s.f)
        return (ka + kb)*cov - (kna + knb)*(1.-cov),ka,kna,kb,knb

rev = ARevPH(1.,100.)
#eta = np.linspace(-1,-0.001,100)
eta = 0.0
#pH = 7.
#pH = 0.
pH = np.linspace(-5,5,100)
r,a,na,b,nb = rev.rate(eta,pH,0.5)
#print r,a,na,b,nb

#log_r = np.log10(abs(r))
#plt.plot(pH,log_r)
plt.plot(pH,r)
plt.show()
