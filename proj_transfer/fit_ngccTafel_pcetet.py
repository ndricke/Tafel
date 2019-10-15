"""
Fits NGCC Tafel data with a pcet2 mechanism

TODO:
Make necessary changes to RevPcetEtScale
Incorporate new mechanism into rest of script
"""

import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import rcParams
import scipy.stats as scst
import scipy.optimize as opt

font = {'size':18}
mpl.rc('font',**font)
rcParams['figure.figsize'] = 14,14


def pcetEt(V, pH, k1o, k2o, kn1o, kn2o, include_reverse=True):
    """ 
    Explicit definition of equation for 2 consecutive PCETs
    Parameters:
    V (numpy 1d array or float): voltage over a domain or at a specific value
    include_reverse (bool): neglect reverse reaction. Sometimes important to disable for numerical stability/testing
    the other inputs are rate constants with pH incorporated. k# represents forward step, kn# reps. backward
    """

    f = 38.949
    expf = np.exp(-0.5*f*V)
    expb = 1./expf
    k1, k2 = k1o*expf, k2o*expf
    kn1, kn2 = kn1o*expb, kn2o*expb
    if include_reverse:
        return (k1*k2 - kn1*kn2)/(k1 + k2 + kn1 + kn2)
    else:
        return (k1*k2)/(k1 + k2 + kn1 + kn2)




n_dom = 100 # number of points for plotting the fitted function

opt_func = pcet2explicit
p_guess = np.array([1.28310184, 1.4140497, 0.43930788, 1.05681383])
p_bounds = ((0., 0., 0., 0.),(np.inf, np.inf, np.inf, np.inf))

## the image for this data is in the same folder as the read_csv
## x: log(J), y: E/V vs NHE
data = pd.read_csv('data_csv/NGCC_Tafel_data.csv', names=['logJ','VvsRHE'])

V_data = np.array(data["VvsRHE"])
J_data = 10.**np.array(data["logJ"]) # shouldn't go to negative infinity during a sign change

# Create Voltage (column 1) and pH (column 2) ranges for plotting the function over a continuous range
V_range = np.linspace(0.6, 0.9, n_dom)

popt, pcov = opt.curve_fit(opt_func, V_data, J_data, p0=p_guess, bounds=p_bounds)

mech_rate = opt_func(V_range, *popt)
print(mech_rate.shape)

## Plot tafel fit and data
plt.scatter(V_range, np.log10(np.abs(mech_rate)), label='Fit: RevPcet2')
plt.scatter(V_data, data['logJ'])

plt.legend()
plt.xlabel("V vs SHE")
plt.ylabel("Log(J)")

plt.show()





