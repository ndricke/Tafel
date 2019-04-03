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

font = {'size':18}
mpl.rc('font',**font)
rcParams['figure.figsize'] = 14,14

n_dom = 100 # number of points for plotting the fitted function
mech = sys.argv[1]

if mech == "pcet2":
    opt_mech = FGH.RevPcet2Scale(dG=-2.46)
    p_guess = (-1.23, 0.01, 0.01, 3*10.**6)
    p_bounds = ((-2.5, 0., 0., 0.),(-0.0000001,np.inf,np.inf,np.inf))

elif mech == "pcetet":
    opt_mech = FGH.RevPcetEtScale(dG=-2.46)
    p_guess = (-1.23, 0.01, 0.01, 3*10.**6)
    p_bounds = ((-2.5, 0., 0., 0.),(-0.0000001,np.inf,np.inf,np.inf))

## the image for this data is in the same folder as the read_csv
## x: log(J), y: E/V vs NHE
data = pd.read_csv('~/work/tafel/tafel_data/NGCC_Tafel_data.csv', names=['logJ','VvsRHE'])
pH_data = pd.read_csv("~/work/tafel/tafel_data/NGCC_pHvsOnset_data.csv", names=["pH", "OnsetV"])

## Convert V vs pH data to RHE, and append it to the tafel data
pH_len = len(pH_data["pH"])
x_data_pH = np.zeros((2,pH_len))
x_data_pH[0,:] = pH_data["OnsetV"] + 0.059*pH_data["pH"] # convert SHE to RHE
x_data_pH[1,:] = pH_data["pH"]

y_data_pH = np.ones(pH_len)*10**-6

x_data = np.array([data["VvsRHE"],data["VvsRHE"]*0.]) #set second column to 0 for pH
y_data = np.array(data["logJ"]) # the data on the plot is actually at -E/V
x_range = np.zeros((2,n_dom))
x_range[0,:] = np.linspace(x_data[0,0], x_data[0,-1], n_dom)

print(x_data)
print(y_data)

popt, pcov = opt.curve_fit(opt_mech.func, x_data, y_data, p0=p_guess, bounds=p_bounds)
mech_rate = opt_mech.func(x_range, *popt)

print(popt) #, pcov)

#fig, (ax0,ax1) = plt.subplots(nrows=2)
fig, ax0 = plt.subplots(nrows=2)

## Plot tafel fit and data
ax0.plot(mech_rate, x_range[0,:], label='Fit: RevPcet2')
ax0.scatter(data['logJ'], data['VvsRHE'])
ax0.legend()
ax0.ylabel("V vs RHE")
ax0.xlabel("Log(J)")

## Plot V vs pH data
#ax1.plot(mech_rate, x_range[0,:], label='Fit: RevPcet2')
#ax1.scatter(data['logJ'], data['VvsRHE'])
#ax1.legend()
#ax1.ylabel("V vs RHE")
#ax1.xlabel("Log(J)")

plt.show()





#plt.savefig("%s_pca%s.png" % (catalyst, pca_n), transparent=True, bbox_inches='tight', pad_inches=0.05)
