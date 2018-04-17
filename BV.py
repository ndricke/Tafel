import scipy as sc
import scipy.optimize as opt
import numpy as np

"""
The simplest electrochemical mechanism, the Butler-Volmer reversible 1-step electron transfer with pH dependence
This will be used for testing simTafel, as well as for subsequently ensuring we can fit the data for the simplest case
"""
class BV(object):
    def __init__(s, pH=7., a=100., b=1., ep=10.**-4):
        s.f = 38.684
        s.pH = pH
        s.ep = ep #what is the smallest point where they can observe the current?
        s.a = a
        s.b = b

    def rate(s, V): #this hopefully allows the user to scan over pH (or f if temperature were interesting?)
        return s.a*10**(-s.pH)*np.exp(-0.5*s.f*V) - s.b*np.exp(0.5*s.f*V)
#        return H*np.exp(-0.5*f*V)# - np.exp(0.5*f*V)

    def fitRate(s, V, *p):
        s.a, s.b = p
        return s.rate(V)

    #Function to generate a voltage-current tafel plot
    def onsetV(s, V):
        return s.rate(V) - s.ep

    def fitOnset(s, pH_linspace, *p):
        s.a, s.b = p
        
#        s.pH = pH
#        V_root = opt.root(s.onsetV, 1.0, method='hybr')
#        return V_root.x[0]

        onsetV_list = [0.] #start with a guess_V within final list so the code can cleanly use previous solution as initial guess

        for i in pH_linspace:
            s.pH = i
            V_root = opt.root(s.onsetV, onsetV_list[-1], method='hybr')
            onsetV_list.append(V_root.x[0])
        onsetV_list.pop(0) #remove the first element, 0, which was only an initial guess
        return onsetV_list

class Mechanism(object):
"""
concentrations: dictionary of the concentrations of non-intermediate species (pH, [O2], ect)
rate_constants: dictionary of the rate constants for every step in the reaction network
ep: for finding pH/concentration dependence, we must assume an observed constant current with varying concentration
"""
    def __init__(s, concentrations, rate_constants, ep):


##Butler-Volmer, with reversible mechanism
class MultiBV(object):
    def __init__(s, a=0.5, ep=10.**-4, params=(1,1,1)):
        s.f = 38.684
        s.a = 0.5
        s.A, s.B, s.C = params
        s.ep = ep #what is the smallest point where they can observe the current?

    def rate(s, V): 
        fwd = (s.A*10.**(-s.pH) + s.B)*np.exp(-s.a*V*s.f)
        bck = (s.A + s.B*10.**(s.pH))*np.exp((1.-s.a)*V*s.f)*s.C
        return fwd - bck

    def fitRate(s, V, *p):
        s.A, s.B, s.C = p
        return s.rate(V)

#Function to generate a voltage-current tafel plot
    def onsetV(s, V):
        return s.rate(V) - s.ep

    def fitOnset(s, pH_linspace, *p):
        s.A, s.B, s.C = p

        onsetV_list = [0.] #start with a guess_V within final list so the code can cleanly use previous solution as initial guess

        for i in pH_linspace:
            s.pH = i
            V_root = opt.root(s.onsetV, onsetV_list[-1], method='hybr')
            onsetV_list.append(V_root.x[0])
        onsetV_list.pop(0) #remove the first element, 0, which was only an initial guess
        return onsetV_list
