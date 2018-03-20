#1. For the ORR mechanism(s) we currently have, produce continuous lines relating: 
#	pH, O2, voltage as expt variables
#	TOF and Tafel as observables (and cyclic voltammetry?)
#	rate constants etc. as parameters
#i.e. sweep over reasonable values of the expt variables/parameters
#2. Now take one example that looks most like the expt system we have, and take a discrete subset of points, add noise, relate noise to paratmer retrieval fidelity.
#3. At this stage: can we tell from this data-to-parameter analysis the difference between two competing mechanisms?
#4. What is the difference between the "Pchem way" and the way they did this in the sci rep paper? In what range of parameter space would these two analyses produce different outcomes? 

import numpy as np
import matplotlib.pyplot as plt

#f: electrothermochemical constant
#a: symmetry coefficient
#k: rate constants; formalism that k_# is a reverse of k#
#H: [H+]
#O2: P_O2

#Function set for inner sphere mechanism of N-doped graphene ORR systems
def evalInSph(eta1,eta2,H,O2,k1,k_1,k2):
    A = inSphA(eta1,H,k1)
    B = inSphB(eta1,k_1)
    C = inSphC(eta2,k2,O2)
    return inSphere(A,B,C)

def inSphere(A,B,C):
    return A*C/(A+B+C)

#need to change the value of f
def inSphA(eta,H,k,a=0.5,f=1.0):
    return k*H*np.exp(-a*f*eta)

def inSphB(eta,k,a=0.5,f=1.0):
    return k*np.exp((1-a)*f*eta)

def inSphC(eta,k,O2,a=0.5,f=1.0):
    return k*O2*np.exp(-a*f*eta)

#the goal is to evaluate a function like inSphere with a variety of functions plugged in for ABC
#This function plugs specific functions taken as arguments into a general functional
#I think we might even want it to return a function that we can repeatedly evaulate?
#XXX This process is still under consideration; I'm not sure if it will be necessary
def eChemEval(XXX,functional,**kwargs):
    result = functional(kwargs)
    return functional(kwargs)






