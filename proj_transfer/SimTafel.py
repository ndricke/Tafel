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


"""
Simulates tafel data for a single electrochemical mechanism, specified in __init__
The parameters of this mechanism should be set elsewhere, then fed to this class
This class is specifically intended to contain all methods that would actually solve kinetic data,
whereas the actual kinetic formulas should be in classes that inherit from mechanism
Note: consider whether it would be more convenient to split this class between methods that simulate and fit data
"""
class SimTafel(object):
    
    def __init__(s, elec_mech, ep=10.0**-1):
        s.elec_mech = elec_mech #class that contains mechanism + rate equation for a proposed mechanistic graph
        s.ep = ep #value for chosen measure of onset potential. Important for calculating onset potential as a function of pH


#Function to generate a voltage-current tafel plot
    def onsetV(s, V):
        return s.elec_mech.rate(V) - s.ep

    def fitLograte(s, V_linspace, *p0):
        s.elec_mech.setConsts(p0)
        return s.elec_mech.lograte(V_linspace)

    def fitOnsetPH(s, pH_linspace, *p0):
        s.elec_mech.setConsts(p0)
        return s.onsetScanPH(pH_linspace)

#Function to calculate the onset potential, V, as a function of pH 
    def onsetScanPH(s, pH_linspace, guess_V=0.):
        onsetV_list = [guess_V] #start with a guess_V within final list so the code can cleanly use previous solution as initial guess

        for i in pH_linspace:
            s.elec_mech.pH = i
            V_root = opt.root(s.onsetV, onsetV_list[-1], method='hybr')
            onsetV_list.append(V_root.x[0])
        onsetV_list.pop(0) #remove the first element, 0, which was only an initial guess
        #print(onsetV_list)
        return onsetV_list

    def findOnsetV(s, onset_J, onsetV_guess=0.):
        """Find the onset potential for a given current
        input: onset_J (float) -- target onset current
        returns: V (float) -- potential that yields onset_J
        """
        s.ep = onset_J
        V_root = opt.root(s.onsetV, onsetV_guess, method='hybr')
        return V_root


    def onsetGradPH(s, pH_linspace, guess_V=1.029):
        onsetV_list = s.onsetScanPH(pH_linspace, guess_V=0)
        dx = pH_linspace[1] - pH_linspace[0]
        gradV = np.gradient(onsetV_list, dx) #the pH[1]-pH[0] is the pH spacing. nernstian d(V_onset)/d(pH) is 59mV/pH
        return onsetV_list, gradV

    def plotTafel(s, V_domain, input_ylim=[-0.2,0.2]):
        I = s.elec_mech.rate(V_domain)
        logI = np.log10(np.abs(I))
        dlogIdV = np.gradient(logI, V_domain[1]-V_domain[0]) #inverse of Tafel slope; calc this way bc V_dom indep var w/ even spacing
        s.dVdlogI = 1./dlogIdV 

        fig, (ax0,ax1) = plt.subplots(nrows=2)
        ax0.plot(V_domain, logI)
        ax1.plot(V_domain, s.dVdlogI)
        plt.xlabel('V')
        ax1.set_ylim(input_ylim)
        ax0.set_xlim([V_domain[0], V_domain[-1]])
        ax1.set_xlim([V_domain[0], V_domain[-1]])
        ax0.set_ylabel('log(I)')
        ax1.set_ylabel('$\partial$V/$\partial$log(I)')
        return fig, ax0, ax1

    def plotdVdPH(s, pH_dom, cut=2):
        onset_Vs, grad_Vs = s.onsetGradPH(pH_dom) 
        fig, (ax0,ax1) = plt.subplots(nrows=2)
        ax0.plot(pH_dom[cut:], onset_Vs[cut:])
        ax1.plot(pH_dom[cut:], grad_Vs[cut:])
        ax0.set_xlim([pH_dom[0],pH_dom[-1]])
        ax1.set_xlim([pH_dom[0],pH_dom[-1]])
        plt.xlabel('pH')
        ax0.set_ylabel(r'V')
        ax1.set_ylabel('$\partial$V/$\partial$pH')
        return fig, ax0, ax1


class FitTafel(SimTafel):

    def __init__(s, elec_data):
        s.elec_data = elec_data

    #Fit the variables for a given mechanism to elec_data
    def fitMech(s, elec_mech):
        pass
        



if __name__ == "__main__":
    
    #import BV
    from ElecMech import BV
    bvm = BV((15., 82.), 1.)


#    V = np.linspace(-0.2,0.2,400)
#    r = bv.rate(V)
#    popt, pcov = opt.curve_fit(bv.fitRate, V, r, p0=(1.0,1.0))    
#    print(popt)

##    plt.plot(V,np.log(r))
#    plt.plot(V,r)


    dom = 400
    pH_dom = np.linspace(3,14,dom)
##    simBV = SimTafel(bv, ep=0)
    bvm = BV((15., 82.), 1.)
    simBV = SimTafel(bvm)
    onset_Vs, grad_Vs = simBV.onsetGradPH(pH_dom)

    noise = (np.random.rand(dom)-0.5)*0.05
    onV_noise = onset_Vs + noise

#    popt, pcov = opt.curve_fit(bv.fitOnset, pH_dom, onset_Vs, p0=(1.,1.))    
    popt, pcov = opt.curve_fit(simBV.fitOnsetPH, pH_dom, onV_noise, p0=(1.,1.))    
    print(popt)

    bv2 = BV((popt[0],popt[1]),1.)
    simBV2 = SimTafel(bv2)
    onset_Vs2, grad_Vs2 = simBV2.onsetGradPH(pH_dom)

    

#    d = {'pH':pH_dom, 'Onset Potential':onset_Vs}
#    pH_onsetV = pd.DataFrame(d)
#    print(pH_onsetV)

    plt.plot(pH_dom, onV_noise)
    plt.plot(pH_dom, onset_Vs2)

    plt.show()

















