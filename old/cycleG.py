import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import scipy.optimize as opt

#reversible cycle via 3 steps, using a general formula with appending pH and potential dependencies
class Rcyc3(object):
    f = 38.949
    ep = 10.**-3
    def __init__(s,I,T):
        s.I = I; s.T = T #all energies are defined in units of kT
        s.n = len(I)
        s.den = np.zeros((s.n,s.n))
        s.dG_r = s.I[-1] - s.I[0]
        s.pHdep = np.array([[0,2,1],[1,0,0],[1,1,0]]) #defined for 1st step not pH dependent
        s.Vdep = np.array([[-1.,1.,0.],[-0.5,-0.5,0.5],[0.5,0.5,-0.5]])

    def rate(s,V,H): #V is applied potential
        chem_sum = H**2 #H is actually a stand-in for OH here, for generality later
        num = np.exp(-s.dG_r)*np.exp(-s.f*V) - chem_sum*np.exp(s.f*V)#as we multiply chem_sum by the reverse rate
        for i in range(s.n):
            for j in range(s.n):
                if i-j >= 0: dGij = s.dG_r
                else: dGij = 0
                s.den[i,j] = np.exp(s.T[i]-s.I[j]-dGij) * H**(s.pHdep[i,j]) * np.exp(s.Vdep[i,j]*s.f*V)
        return num/np.sum(s.den)

    def varV(s,V):
        return s.rate(V,s.OH) - s.ep

    def varpH(s,OH):
        return rate(s.V,OH) - s.ep

#Generic reversible cycle, formula derived using steady state approximation, no branching, linear in intermediates
class RcycG(object): 
    f = 38.949 
    ep = 10.**-1
    def __init__(s,ktab,Vdep,pHdep):
        s.k = ktab
        s.n = s.k.shape[1]
        s.Vdep = Vdep
        s.pHdep = pHdep

    def rate(s,V,H): #V is applied potential
        kmod = s.depVPH(V,H)
        num = np.prod(kmod[0,:])-np.prod(kmod[1,:])
        den = np.ones((s.n,s.n))
        for i in range(s.n):
            for j in range(s.n):
                k_list,indx_arr = s.denTerms(i,j)
                for te in range(s.n-1):
                    den[i,j] *= kmod[indx_arr[te],k_list[te]]
        return num/np.sum(den)

    def depVPH(s,V,H):
        kmod = np.zeros((2,s.n))
        OH = 10**-14/H
        kmod[0,:] = s.k[0,:] * H**(s.pHdep[0,:]) * np.exp(-0.5*s.Vdep*s.f*V)
        kmod[1,:] = s.k[1,:] * OH**(s.pHdep[1,:]) * np.exp(0.5*s.Vdep*s.f*V)
        return kmod

    def denTerms(s,i,j): #i defines which rate constant is missing; j defines which intermediate is pointed
        k_list = range(s.n) #list from 0 to n-1 of n terms
        k_list.remove(i) #remove the missing rate constant

        sign_arr = np.ones(s.n-1)#list of indices that determine sign of rate constant. 1 is fwd, -1 is back
        i_j = np.sign(i-j)
        sign_arr[i:j] *= -1 # switch all terms in between i and j
        sign_arr *= np.sign(i-j) # switch all terms if i < j
        indx_arr = [int(c) for c in -0.5*sign_arr+0.5] #convert to 0 for fwd, 1 for back

        return k_list, indx_arr

    def varV(s,V):
        return s.rate(V,s.H) - s.ep

    def varpH(s,H):
        return rate(s.V,H) - s.ep

if __name__ == "__main__":

    pH_range = [11.,13.5]
    h = np.linspace(pH_range[0],pH_range[1],200)
    OH = 10.**(h-14.)

    V_range = np.linspace(-3.,3.,200)

    I = np.array([-1.,1.,-0.5]) 
    ts = np.array([3.,2.,2.]) #transition state energy above previous intermediate
    T = I+ts #energy of TS relative to initial intermediate

    cyc = Rcyc3(I,T)
    #print cyc.rate(0,-1.08)
    #print cyc.den

    #rate_scan = []
    #for i in OH:
    #    rate_scan.append(cyc.rate(0,i))
    #plt.plot(h,rate_scan)
    #plt.show()

    #rate_scan = []
    #OH_const = 10**-3
    #for i in V_range:
    #    rate_scan.append(cyc.rate(i,OH_const))
    #plt.plot(V_range,rate_scan)
    #plt.show()

    #cyc.OH = 10**-3
    #sol = opt.root(cyc.varV,0,method='hybr')
    #print sol

    sol_list = [0] #start with a 0 so the code can cleanly use this as an initial guess
    for i in OH:
        cyc.OH = i
        sol = opt.root(cyc.varV,sol_list[-1],method='hybr')
        sol_list.append(sol.x[0])
    sol_list.pop(0)

    plt.plot(h,sol_list)
    plt.xlabel("pH")
    plt.ylabel("Potential (V)") #when current is observable
    plt.xlim(pH_range)
    plt.show()








