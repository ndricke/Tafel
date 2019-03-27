import sys

import numpy as np
import pandas as pd
#import tafel.DisorderMech as DM
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as scst
import scipy.optimize as opt

import tafel.ElecMech
import fit_goldHER


"""
Script for fitting NGCC tafel and pH data
There is a weak pH dependence, ~0.8 order O2 dependence, and tafel slope shifting from 60-120 mV/dec
The RDS is probably an O2-coupled electron transfer, but there must be a pH dependent step that is likely
mixed in. This could be either a PCET or proton transfer, before or afterward (although I think these are equivalent)
"""

if __name__ == "__main__":
    font = {'size':18}
    mpl.rc('font',**font)

    ## the image for this data is in the same folder as the read_csv
    ## x: log(J), y: E/V vs NHE
    tafel_data = pd.read_csv('~/Pictures/tafel/NGCC_Tafel_data.csv', names=['logJ','VvsRHE'])
    pH_data = pd.read_csv('~/Pictures/tafel/NGCC_pHvsOnset_data.csv', names=['pH','VvsNHE'])
    O2_data = pd.read_csv('~/Pictures/tafel/NGCC_O2vsRate_data.csv', names=['pO2atm','logJ'])

    # check signs for these to make sure they match equations
    y_data = np.array(tafel_data['logJ'])
    n_data = len(y_data)
    pH = 0
    H_conc = 10.**(-pH)
    x_data = np.array([tafel_data['VvsRHE'],H_conc*np.ones(n_data)])
    print(x_data.shape)

    ## We'll use disorder when we're ready
    #sig = 0.1
    #cyc_count = 5
    #deltas, weights = np.polynomial.hermite.hermgauss(cyc_count)  #gaussian quadrature, x=deltas, y=weights
    #weights *= 1./(np.pi**0.5) # the weights from gaussian quadrature actually sum to 2
    #deltas *= 0.5*sig*2**0.5

    # Fit Tafel data with PcetEt
    opt_mech = fit_goldHER.RevPcetEt()
    popt = (-.6,0.0000001,0.000001,-1.23)
    opt_mech.dG = popt[-1]
    opt_mech.genConsts(popt[:-1])
    print(opt_mech.getConsts())



    ##popt, pcov = opt.curve_fit(opt_mech.funcFullE , x_data, y_data, p0=(-0.2,0.001,0.001,0.0), bounds=((-5.0,0.,0.,-5.0),(5.0,np.inf,np.inf,5.0)))
    #mech_rate = opt_mech.funcFullE(x_data, popt[0], popt[1], popt[2], popt[3])

    #opt_mech = fit_goldHER.RevPcetEt(weights=weights, ddGs=deltas, mech="acid")
    #popt, pcov = opt.curve_fit(opt_mech.funcDisFullE, x_data, y_data, p0=(-0.2,0.01,0.01,1.2,0.0), bounds=((-5.0,0.,0.,0.,-5.0),(5.0,np.inf,np.inf,5.0,5.0)))
    #mech_rate = opt_mech.funcDisFullE(x_data, popt[0], popt[1], popt[2], popt[3], popt[4]) #optimizing for disorder makes it real big

    print(opt_mech.rate(0.,1.))

    #print(x_data[0,:].shape)
    #print(y_data.shape)
    #print(mech_rate.shape)
#
    #print(x_data)
    #print()
    #print(y_data)
    #print()
    #print(mech_rate)

    #plt.scatter(x_data[0,:], y_data)
    #plt.scatter(x_data[0,:], mech_rate)
    #plt.show()



##
