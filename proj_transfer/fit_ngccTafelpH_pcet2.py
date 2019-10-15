"""
Fits NGCC Tafel data with a pcet2 mechanism without ignore pH

TODO:
This fit is not very good. I know the mechanism shouldn't be able to fit it very well, but I think there
is probably some issue with the initial guess or some assumption about the parameters that is keeping it
from being able to fit both pieces of data very well.
"""

import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import rcParams
import scipy.stats as scst 
import scipy.optimize as opt 
import EChemLib

font = {'size':18}
mpl.rc('font',**font)
rcParams['figure.figsize'] = 14,14

def pcet2VpH(VpH, k1o, k2o, kn1o, kn2o, pH_mediated="base", include_reverse=False):
    """
    Same idea as the function pcet2explicit, but VpH is a numpy 2d array,
    where the first column is voltage, and the second is pH
    
    Mechanism: 
    R + e- + H2O <--> R-H + OH- | Base Mediated
    R + e- + H+ <--> R-H        | Acid Mediated

    [H+] = 10**(-pH)
    [H+][OH-] = 10**-14
    pH + pOH = 14
    pOH = 14 - pH
    [OH-] = 10**(-pOH) = 10**(pH - 14)

    k1o (float): first step forward rate constant
    k2o (float) second step forward rate constant
    kn1o (float): first step backward rate constant
    kn2o (float) second step backward rate constant
    pH_mediated (str): acid or base. Dictates possible proton sources for mechanism
    include_reverse (bool): neglect reverse reaction
    """
    f = 38.949
    expf = np.exp(-0.5*f*VpH[:,0])
    expb = 1./expf
    
    if pH_mediated == "acid":
        expf *= 10.**(VpH[:,1])
    elif pH_mediated == "base":
        expb *= 10.**(VpH[:,1] - 14.)

    k1, k2 = k1o*expf, k2o*expf
    kn1, kn2 = kn1o*expb, kn2o*expb
    if include_reverse:
        return (k1*k2 - kn1*kn2)/(k1 + k2 + kn1 + kn2)
    else:
        return (k1*k2)/(k1 + k2 + kn1 + kn2)


n_dom = 100 # number of points for plotting the fitted function
pH = 1.
mediation = "acid"

opt_func = pcet2VpH
p_guess = np.array([1.28310184, 1.4140497, 0.43930788, 1.05681383])
p_bounds = ((0., 0., 0., 0.),(np.inf, np.inf, np.inf, np.inf))

## the image for this data is in the same folder as the read_csv
## x: log(J), y: E/V vs NHE
data = pd.read_csv('data_csv/NGCC_Tafel_data.csv', names=['logJ','VvsRHE'])
data_pH = pd.read_csv("data_csv/NGCC_pHvsOnset_data.csv", names = ["pH", "VvsSHE"])  

# Convert RHE to SHE
data["VvsSHE"] = EChemLib.convertHE(np.array(data["VvsRHE"]), pH=13.0, H_electrode="SHE")



V_data = np.array(data["VvsSHE"])
VpH = np.zeros((V_data.shape[0],2))
VpH[:,0] = V_data
VpH[:,1] = pH

VpH_onset = np.array(data_pH)
VpH_onset[:, [0,1]] = VpH_onset[:, [1,0]]

VpH = np.concatenate((VpH, VpH_onset), axis=0)
print("V,          pH")
print(VpH)

J_data = 10.**np.array(data["logJ"]) # shouldn't go to negative infinity during a sign change
J_pH_data = np.ones(data_pH.shape[0])*0.05  # onset potential measured at const 0.05 mA/cm^-2
J_data = np.concatenate((J_data, J_pH_data))

# Create Voltage (column 1) and pH (column 2) ranges for plotting the function over a continuous range
V_range = np.linspace(0.6, 0.9, n_dom)
VpH_range = np.zeros((n_dom, 2))
VpH_range[:,0] = V_range
VpH_range[:,1] = pH

popt, pcov = opt.curve_fit(opt_func, VpH, J_data, p0=p_guess, bounds=p_bounds)

mech_rate = opt_func(VpH_range, *popt)
mech_rate_constJ = opt_func(VpH_onset, *popt)
print(mech_rate.shape)

fig, (ax0, ax1) = plt.subplots(nrows=2)

## Plot tafel fit and data
ax0.scatter(V_range, np.log10(np.abs(mech_rate)), label='Fit: Pcet2')
ax0.scatter(V_data, data['logJ'])
ax0.legend()

ax1.scatter(data_pH["VvsSHE"], np.log10(np.abs(mech_rate_constJ)), label="Fit: Pcet2")
ax1.scatter(data_pH["VvsSHE"], np.log10(np.abs(J_pH_data)))
ax0.legend()

plt.xlabel("V vs SHE")
plt.ylabel("Log(J)")
#fig.set_size(16., 8., forward=True)

plt.show()





