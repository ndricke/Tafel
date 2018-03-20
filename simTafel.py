"""
03/19/18
By: Nathan Ricke
A class that serves as a library for generating simulated current-voltage data (tafel plots)
The structure is a class so that a single set of initial parameters may be shared among the varying data generating functions

Data needed by each generating function:
1. A voltage scan range + datapoints (np.linspace)
2. The voltage-current equation derived from the mechanism's graph

It would be sensible for a separate piece of software to convert mechanism graphs into voltage-current equations.
Given that the software that generates mechanisms could change form appreciably, it probably makes the most sense to have that
within a separate file that we can import here. 

Now, one of the things we would be interested in fitting is the pH relationships
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import scipy.optimize as opt

import cycleG

"""
The simplest electrochemical mechanism, the Butler-Volmer reversible 1-step electron transfer with pH dependence
This will be used for testing simTafel, as well as for subsequently ensuring we can fit the data for the simplest case
"""
class BV(object):
    def __init__(s, pH=7., fwd_scale=0.01, ep=10.**-4):
        s.f = 38.684
        s.pH = pH
        s.ep = ep #what is the smallest point where they can observe the current?
        s.fwd_scale = fwd_scale

    def rate(s, V): #this hopefully allows the user to scan over pH (or f if temperature were interesting?)
        return s.fwd_scale*10**(-s.pH)*np.exp(-0.5*s.f*V) - np.exp(0.5*s.f*V)
#        return H*np.exp(-0.5*f*V)# - np.exp(0.5*f*V)



"""
Simulates tafel data for a single electrochemical mechanism, specified in __init__
The parameters of this mechanism should be set elsewhere, then fed to this class
"""

class simTafel(object):
    
    def __init__(s, elec_mech, ep=10.0**-4):
        s.elec_mech = elec_mech #class that contains mechanism + rate equation for a proposed mechanistic graph
        s.ep = ep #value for chosen measure of onset potential. Important for calculating onset potential as a function of pH


#Function to generate a voltage-current tafel plot
    def onsetV(s, V):
        return s.elec_mech.rate(V) - s.ep

#Function to calculate the onset potential, V, as a function of pH 
    def OnsetScanPH(s, pH_linspace, guess_V=0.):
        onsetV_list = [guess_V] #start with a guess_V within final list so the code can cleanly use previous solution as initial guess

        for i in pH_linspace:
            s.elec_mech.pH = i
            V_root = opt.root(s.onsetV, onsetV_list[-1], method='hybr')
            onsetV_list.append(V_root.x[0])
        onsetV_list.pop(0) #remove the first element, 0, which was only an initial guess

        gradV = np.gradient(onsetV_list, pH_linspace[1]-pH_linspace[0]) #the pH[1]-pH[0] is the pH spacing. nernstian d(V_onset)/d(pH) is 59mV/pH
        return onsetV_list



if __name__ == "__main__":
    
    bv = BV()


#    V = np.linspace(-0.2,0.2,400)
#    r = bv.rate(V)
#
##    plt.plot(V,np.log(r))
#    plt.plot(V,r)

    pH_dom = np.linspace(3,14,400)
    simBV = simTafel(bv)
    onset_Vs = simBV.OnsetScanPH(pH_dom)

    plt.plot(pH_dom, onset_Vs)

    plt.show()
