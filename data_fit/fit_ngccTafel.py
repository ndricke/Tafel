import sys

import numpy as np
import pandas as pd
#import tafel.DisorderMech as DM
import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import rcParams

import scipy.stats as scst
import scipy.optimize as opt


import tafel.data_fit.fit_goldHER as FGH


def pcet2explicit(V, k1o, k2o, kn1o, kn2o):
    f = 38.949
    expf = np.exp(-0.5*f*V)
    expb = 1./expf
    k1, k2 = k1o*expf, k2o*expf
    kn1, kn2 = kn1o*expb, kn2o*expb
    #return (k1*k2 - kn1*kn2)/(k1 + k2 + kn1 + kn2)
    return (k1*k2)/(k1 + k2 + kn1 + kn2)

def pcet2VpH(VpH, k1o, k2o, kn1o, kn2o):
    f = 38.949
    expf = np.exp(-0.5*f*V[:,0])
    expb = 1./expf
    expfpH = expf*10.**(VpH[:,1] - 14.)
    k1, k2 = k1o*expfpH, k2o*expfVpH
    kn1, kn2 = kn1o*expb, kn2o*expb
    #return (k1*k2 - kn1*kn2)/(k1 + k2 + kn1 + kn2)
    return (k1*k2)/(k1 + k2 + kn1 + kn2)


font = {'size':18}
mpl.rc('font',**font)
rcParams['figure.figsize'] = 14,14

n_dom = 1000 # number of points for plotting the fitted function
mech = sys.argv[1]

if mech == "pcet2":
    opt_mech = FGH.RevPcet2Scale(dG=-2.46)
    #p_guess = (-1.23, 2.57e-02, 4.57e-02, 11.)
    #p_guess = (-1.23, 0.01, 0.01, 3*10.**12)
    #p_guess = (-1.23, 0.01, 0.01, 3*10.**-6)
    #p_bounds = ((-2.5, 0., 0., 0.),(-0.0000001,np.inf,np.inf,np.inf))
elif mech == "ktest":
    opt_func = pcet2explicit
    #p_guess = np.array([1.28310184, 1.4140497, 0.43930788, 0.05681383])
    p_guess = np.array([1.28310184, 1.4140497, 0.43930788, 0.05681383])
    p_guess[3] /= 1000.
    #p_guess = (12.8310184, 14.140497, 0.43930788, 0.05681383)
    p_bounds = ((0., 0., 0., 0.),(np.inf,np.inf,np.inf,np.inf))
elif mech == "ktestpH":
    opt_func = pcet2VpH
    #p_guess = np.array([1.28310184, 1.4140497, 0.43930788, 0.05681383])
    p_guess = np.array([1.28310184, 1.4140497, 0.43930788, 0.05681383])
    p_guess[3] /= 1000.
    #p_guess = (12.8310184, 14.140497, 0.43930788, 0.05681383)
    p_bounds = ((0., 0., 0., 0.),(np.inf,np.inf,np.inf,np.inf))
elif mech == "pcetet":
    opt_mech = FGH.RevPcetEtScale(dG=-2.46)
    p_guess = (-1.23, 0.01, 0.01, 3*10.**6)
    p_bounds = ((-2.5, 0., 0., 0.),(-0.0000001,np.inf,np.inf,np.inf))

## the image for this data is in the same folder as the read_csv
## x: log(J), y: E/V vs NHE
data = pd.read_csv('~/work/tafel/tafel_data/NGCC_Tafel_data.csv', names=['logJ','VvsRHE'])
pH_data = pd.read_csv("~/work/tafel/tafel_data/NGCC_pHvsOnset_data.csv", names=["pH", "OnsetV"])

## Convert V vs pH data is already in SHE
pH_len = len(pH_data["pH"])
x_data_pH = np.zeros((2,pH_len))
x_data_pH[0,:] = pH_data["OnsetV"] # this is already vs NHE
x_data_pH[1,:] = pH_data["pH"]

y_data_pH = np.ones(pH_len)*10**-6

x_data = np.zeros((len(data),2))
x_data[:,0] = data["VvsRHE"]
pH = 13. #the tafel data is collected all at pH 13
x_data[0,:] -= 0.059 * pH # convert Tafel data RHE to NHE for direct comparison
y_data = 10.**np.array(data["logJ"]) # convert y_data to J rather than log10(J)
x_range = np.zeros((2,n_dom))
x_range[0,:] = np.linspace(x_data[0,0]-0.35, x_data[0,-1]+0.35, n_dom)
x_range[1,:] = np.ones(n_dom)*13.

print("x data: ", x_data)
print("y data: ", y_data)

popt, pcov = opt.curve_fit(opt_mech.func, x_data, y_data, p0=p_guess, bounds=p_bounds)
#popt, pcov = opt.curve_fit(opt_func, x_data[0,:].flatten(), y_data, p0=p_guess, bounds=p_bounds)
#popt = p_guess
#mech_rate = opt_func(x_range[0,:].flatten(), *popt)

#popt = p_guess
mech_rate = opt_mech.func(x_range, *popt)

print("Optimized Parameters: ", popt) #, pcov)
print("Parameter covariance matrix: ")
#print(pcov)

#fig, (ax0,ax1) = plt.subplots(nrows=2)
#fig, ax0 = plt.subplots(nrows=2)

## Instanteneous Tafel slope of data
logJ = data['logJ'].values
mV = x_data[0,:]
print("mV")
print(mV)
print("logJ")
print(logJ)
slope = (mV[1:] - mV[:-1])/(logJ[1:] - logJ[:-1])
print("mV/logJ")
print(slope)

## Plot tafel fit and data
#plt.plot(np.log10(np.abs(mech_rate)), x_range[0,:], label='Fit: RevPcet2')
plt.plot(x_range[0,:], np.log10(np.abs(mech_rate)), label='Fit: RevPcet2')
#plt.scatter(data['logJ'], data['VvsRHE'])
plt.scatter(x_data[0,:], data['logJ'])
plt.legend()
plt.xlabel("V vs SHE")
plt.ylabel("Log(J)")
#plt.xlim(np.min(data['logJ'])-0.1, np.max(data['logJ']+0.1))
#plt.ylim(np.min(x_data[1,:])-0.1, np.max(x_data[0,:])+0.1)

## Plot V vs pH data
#ax1.plot(mech_rate, x_range[0,:], label='Fit: RevPcet2')
#ax1.scatter(data['logJ'], data['VvsRHE'])
#ax1.legend()
#ax1.ylabel("V vs RHE")
#ax1.xlabel("Log(J)")

plt.show()





#plt.savefig("irpcet2_gccTafel.png", transparent=True, bbox_inches='tight', pad_inches=0.05)
#plt.savefig("%s_pca%s.png" % (catalyst, pca_n), transparent=True, bbox_inches='tight', pad_inches=0.05)
