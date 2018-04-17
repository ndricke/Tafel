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

import sys

import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import scipy.optimize as opt
import pandas as pd

import cycleG
from BV import BV, MultiBV




"""
Simulates tafel data for a single electrochemical mechanism, specified in __init__
The parameters of this mechanism should be set elsewhere, then fed to this class
"""
class SimTafel(object):
    
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
        return onsetV_list, gradV


class FitTafel(object):

    def __init__(s, elec_data):
        s.elec_data = elec_data

    #Fit the variables for a given mechanism to elec_data
    def fitMech(s, elec_mech):
        pass
        



if __name__ == "__main__":
    
    bv = BV(a=15,b=82.2)


#    V = np.linspace(-0.2,0.2,400)
#    r = bv.rate(V)
#    popt, pcov = opt.curve_fit(bv.fitRate, V, r, p0=(1.0,1.0))    
#    print(popt)

##    plt.plot(V,np.log(r))
#    plt.plot(V,r)


    dom = 400
    pH_dom = np.linspace(3,14,dom)
##    simBV = SimTafel(bv, ep=0)
    simBV = SimTafel(bv)
    onset_Vs, grad_Vs = simBV.OnsetScanPH(pH_dom)

    noise = (np.random.rand(dom)-0.5)*0.05
    onV_noise = onset_Vs + noise

#    popt, pcov = opt.curve_fit(bv.fitOnset, pH_dom, onset_Vs, p0=(1.,1.))    
    popt, pcov = opt.curve_fit(bv.fitOnset, pH_dom, onV_noise, p0=(1.,1.))    
    print(popt)

    bv2 = BV(a=popt[0], b=popt[1])
    simBV2 = SimTafel(bv2)
    onset_Vs2, grad_Vs2 = simBV2.OnsetScanPH(pH_dom)

    

#    d = {'pH':pH_dom, 'Onset Potential':onset_Vs}
#    pH_onsetV = pd.DataFrame(d)
#    print(pH_onsetV)

    plt.plot(pH_dom, onV_noise)
    plt.plot(pH_dom, onset_Vs2)

    plt.show()

















