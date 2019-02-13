import sys

import numpy as np
import pandas as pd
#import tafel.DisorderMech as DM
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as scst
import scipy.optimize as opt

import tafel.ElecMech

font = {'size':18}
mpl.rc('font',**font)

## the image for this data is in the same folder as the read_csv
## x: log(J), y: E/V vs NHE
data = pd.read_csv('~/pymod/tafel/data/Jackson_HER_data.csv', names=['logJ','VvsNHE'])

"""
this data actually has a 3rd column, which is pH. Need to add that data to this dataframe
every pH unit has 6 data points, except pH = 3
the data is sorted in order of ascending pH
the data also skips pH = 5. Overall there are 65 data points
"""

class RevPcet2(tafel.ElecMech.ElecMech):

    def __init__(self, weights=None, ddGs=None, mech="acid", dG=0., a=0.5):
        """
        Parameters:
        weights: (list) population of each cycle (normalized to 1). fraction of active sites for each mech
        dG: (float) free energy of overall cycle
        """

        self.kab = 1. #this is to be set explicitly if acid and base mechanisms have different barriers
        self.kab2 = 1.
        super().__init__(dG=dG, a=a)
        if mech == "base" or mech == "both":
            self.mech = mech
        elif mech == "acid":
            self.mech = "acid"
        else:
            raise ValueError("%s is not a valid PCET mechanism" % (str(mech)) )

        if ddGs is None and weights is None:
            self.ddGs = [0.]
            self.weights = [1.]
        else:
            self.weights = weights
            self.ddGs = ddGs

    def rate(self, V, pH):
        H = 10.**(-pH)
        #k1, k2 = self.k1*H*np.exp(-0.5*self.f*V), self.k2*H*np.exp(-0.5*self.f*V)
        #kn1, kn2 = self.kn1*np.exp(0.5*self.f*V), self.kn2*np.exp(0.5*self.f*V)
        k1, kn1 = self.PCET((self.k1,self.kn1, self.kab), V=V, H=H, mech=self.mech)
        k2, kn2 = self.PCET((self.k2,self.kn2, self.kab2), V=V, H=H, mech=self.mech)
        return (k1*k2 - kn1*kn2)/(k1 + k2 + kn1 + kn2)

    def setConsts(self, p):
        self.k1, self.k2, self.kn1, self.kn2 = p

    def setConc(self, p):
        self.pH = p

    """Generate the rate constants for this reaction from intermediate and TS energies"""
    def genConsts(self, dGi_Ti):
        self.k1, self.kn1 = self.gen_k_edge(dGi_Ti[0], dGi_Ti[1])
        self.k2, self.kn2 = self.gen_k_edge(self.dG-dGi_Ti[0], dGi_Ti[2]) #for HER, dG of the cycle is 0 vs NHE

    def func(self, VpH, dG1, T1, T2):
        self.genConsts((dG1, T1, T2))
        rate = self.rate(VpH[0,:], VpH[1,:])
        return np.log10(rate)

    def disorderFunc(self, VpH, dG1, T1, T2, sig):
        rate = 0.
        mod_ddGs = self.ddGs*0.5*sig*2**0.5
        for i in range(len(self.weights)):
            self.genConsts((self.dG, dG1+mod_ddGs[i], T1, T2))
            rate += self.weights[i]*self.rate(VpH[0,:], VpH[1,:])
        return np.log10(rate)

    def funcDisAB(self, VpH, dG1, T1, T2, sig, kab1, kab2):
        self.kab = kab1
        self.kab2 = kab2
        return self.disorderFunc(VpH, dG1, T1, T2, sig)

    def funcDisFullE(self, VpH, dG1, T1, T2, sig, dG):
        self.dG = dG
        return self.disorderFunc(VpH, dG1, T1, T2, sig)

    def funcDisFullEAB(self, VpH, dG1, T1, T2, sig, dG, kab):
        self.dG = dG
        self.kab = kab
        return self.disorderFunc(VpH, dG1, T1, T2, sig)

class RevPcetEt(RevPcet2):

    def rate(self, V, pH):
        H = 10.**(-pH)
        k1, kn1 = self.PCET((self.k1, self.kn1, self.kab), V=V, H=H, mech=self.mech)
        k2, kn2 = self.k2*np.exp(-0.5*self.f*V), self.kn2*np.exp(0.5*self.f*V)
        return (k1*k2 - kn1*kn2)/(k1 + k2 + kn1 + kn2)




pH_values = [1,2,3,4,6,7,8,9,10,11,12] #1-12 except 5
pH_list = []
for i in pH_values:
    for j in range(6):
        pH_list.append(i)
pH_list.remove(3) #pH=3 only has 5 data points

data = data.assign(pH=pH_list)
#print(data[:29])
#print(data[29:])
data = data[data["pH"] != 4]

fig, ax = plt.subplots()

x_data = np.array([-1.*data["VvsNHE"],data["pH"]])
y_data = np.array(data["logJ"]) # the data on the plot is actually at -E/V

sig = 0.6
cyc_count = 7
deltas, weights = np.polynomial.hermite.hermgauss(cyc_count)  #gaussian quadrature, x=deltas, y=weights
weights *= 1./(np.pi**0.5) # the weights from gaussian quadrature actually sum to 2
deltas *= 0.5*sig*2**0.5

### Plot of MJdata with 2 reversible PCET's (no disorder)
#opt_mech = RevPcet2()
#popt, pcov = opt.curve_fit(opt_mech.func, x_data, y_data, p0=(-0.2,0.01,0.01), bounds=((-2.,0.,0.),(-0.0000001,np.inf,np.inf)))
#print(popt, pcov)
#mech_rate = opt_mech.func(x_data, popt[0], popt[1], popt[2])
#plt.scatter(mech_rate, -1.*x_data[0,:], label='Fit: no disorder')

### Plot of MJdata with 2 reversible PCET's (disorder on first step)
#opt_mech = RevPcet2(weights=weights, ddGs=deltas)
#popt, pcov = opt.curve_fit(opt_mech.disorderFunc, x_data, y_data, p0=(-0.2,0.01,0.01), bounds=((-np.inf,0.,0.),(np.inf,np.inf,np.inf)))
#mech_rate = opt_mech.disorderFunc(x_data, popt[0], popt[1], popt[2])
#plt.scatter(mech_rate, -1.*x_data[0,:], label='Fit: Disorder')
#mech_rate_mf = opt_mech.func(x_data, popt[0], popt[1], popt[2])
#plt.scatter(mech_rate_mf, -1.*x_data[0,:], label='NoDis w/ dis consts)')

### Plot MJdata fit with PCET-ET (no disorder) --> the pH dependence on Pcet2 is too high, so this curbs it
#opt_mech = RevPcetEt()
#popt, pcov = opt.curve_fit(opt_mech.func, x_data, y_data, p0=(-0.2,0.01,0.01), bounds=((-2.,0.,0.),(-0.0000001,np.inf,np.inf)))
#mech_rate = opt_mech.func(x_data, popt[0], popt[1], popt[2])


#opt_mech = RevPcetEt(weights=weights, ddGs=deltas)
### Plot MJdata fit with PCET-ET (disorder)
#popt, pcov = opt.curve_fit(opt_mech.disorderFunc, x_data, y_data, p0=(-0.2,0.01,0.01,0.4), bounds=((-5.0,0.,0.,0.),(5.0,np.inf,np.inf,5.0)))
##mech_rate = opt_mech.disorderFunc(x_data, popt[0], popt[1], popt[2])
#mech_rate = opt_mech.disorderFunc(x_data, popt[0], popt[1], popt[2], popt[3]) #optimizing for disorder makes it real big


### Plot MJdata with total energy of mechanism as a parameter
#popt, pcov = opt.curve_fit(opt_mech.funcDisFullE, x_data, y_data, p0=(-0.2,0.01,0.01,1.2,0.0), bounds=((-5.0,0.,0.,0.,-5.0),(5.0,np.inf,np.inf,5.0,5.0)))
#mech_rate = opt_mech.funcDisFullE(x_data, popt[0], popt[1], popt[2], popt[3], popt[4]) #optimizing for disorder makes it real big

### Plot MJdata with acid-base mechanism switch
#opt_mech = RevPcetEt(weights=weights, ddGs=deltas, mech="both", dG=0.)
#popt, pcov = opt.curve_fit(opt_mech.funcDisFullEAB, x_data, y_data, p0=(-0.2,0.01,0.01,0.4,0.0,5.0), \
#                           bounds=((-5.0,0.,0.,0.,-5.0,0.),(5.0,np.inf,np.inf,5.0,5.0,np.inf)))
#mech_rate = opt_mech.funcDisFullEAB(x_data, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5]) #optimizing for disorder makes it real big


### Plot MJdata with acid-base mechanism switch on Pcet2
opt_mech = RevPcet2(weights=weights, ddGs=deltas, mech="both", dG=-0.01)
popt, pcov = opt.curve_fit(opt_mech.funcDisFullEAB, x_data, y_data, p0=(-0.1,0.01,0.01,0.2,1.0,1.0), \
                           bounds=((-5.0,0.,0.,0.,0.,0.),(5.0,np.inf,np.inf,5.0,np.inf,np.inf)))
mech_rate = opt_mech.funcDisFullEAB(x_data, popt[0], popt[1], popt[2], popt[3], popt[4], popt[5]) #optimizing for disorder makes it real big

print(popt) #, pcov)
plt.scatter(mech_rate, -1.*x_data[0,:], label='Fit: Disorder')


### Plot experimental MJdata
slope_list, intercept_list, r_val_list = [], [], []
x = np.array([-4.30, -2.6])
for pH_val in pH_values:
    plt.scatter(data[data['pH']==pH_val]['logJ'], data[data['pH']==pH_val]['VvsNHE'], label=pH_val)

    ### Fit a line to each pH value, plot it, and save relevant values
    #slope, intercept, r_val, p_val, std_err = scst.linregress(data[data['pH']==pH_val]['logJ'], data[data['pH']==pH_val]['VvsNHE'])
    ##plt.plot(x, slope*x + intercept)
    #slope_list.append(slope)
    #intercept_list.append(intercept)
    #r_val_list.append(r_val)


## Return slopes of lines
#print(slope_list)
#print(intercept_list)
#print(r_val_list)

fig.set_size_inches(11.,11., forward=True)



## Plot of MJdata with linear regression fits
plt.ylabel("-E/V vs NHE")
plt.xlabel(r"log(J/Acm$^{-2}$)")
#plt.legend(loc=3)
#plt.savefig("MJdata_linreg.png", transparent=True, bbox_inches='tight', pad_inches=0.02)
#plt.savefig("MJdata_revpcet2.png", transparent=True, bbox_inches='tight', pad_inches=0.02)
#plt.savefig("MJdata_revpcet2_disorder.png", transparent=True, bbox_inches='tight', pad_inches=0.02)
#plt.savefig("MJdata_RevPcetEt.png", transparent=True, bbox_inches='tight', pad_inches=0.02)
#plt.savefig("MJdata_RevPcetEt_disorder.png", transparent=True, bbox_inches='tight', pad_inches=0.02)
#plt.savefig("MJdata_RevPcetEt_FreeDisorder.png", transparent=True, bbox_inches='tight', pad_inches=0.02)
#plt.savefig("MJdata_RevPcetEt_FreeDisorder_dGfixed.png", transparent=True, bbox_inches='tight', pad_inches=0.02)
#plt.savefig("MJdata_RevPcetEt_FreeDisorder_High.png", transparent=True, bbox_inches='tight', pad_inches=0.02)
#plt.savefig("MJdata_RevPcetEt_FreeDisorder_AB.png", transparent=True, bbox_inches='tight', pad_inches=0.02)

## Plot of linear regression slopes as a function of pH
#plt.scatter(pH_values, slope_list, s=100)
#plt.plot(pH_values, slope_list)
#plt.xlabel("pH")
#plt.ylabel("Linear Regression Slope")
#plt.savefig("MJdata_linreg_slopes.png", transparent=True, bbox_inches='tight', pad_inches=0.02)

plt.show()
