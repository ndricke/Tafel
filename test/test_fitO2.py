import sys

import scipy as sc
from scipy import stats
import scipy.optimize as opt
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score

import ElecMech
import intrinsic_het as IH
import SimTafel

matplotlib.rcParams.update({'font.size': 22})

class DataFit:
    def __init__(s, emech, V, pO2dom, Jdata, disorder=False):
        s.emech = emech
        s.Jdata = Jdata
        s.pO2dom = pO2dom
        s.V = V
        s.disorder = disorder

    def evalfit(s, kconsts):
        rate_arr = s.O2vsRate(kconsts)
        return np.linalg.norm(rate_arr - s.Jdata)

    def O2vsRate(s, kconsts, pO2 = None):
        if pO2 is None:
            pO2 = s.pO2dom
        if s.disorder:
            ratefunc = s.emech.rate
        else:
            ratefunc = s.emech.baseRate
        s.emech.setConsts(kconsts)
        rate_arr = np.zeros(len(pO2))
        for i, point in enumerate(pO2):
            s.emech.O2 = pO2[i]
            rate_arr[i] = ratefunc(s.V)
        return rate_arr

    def evalfitDis(s, kconsts):
        rate_arr = s.O2vsRate(kconsts)
        return np.linalg.norm(rate_arr - s.Jdata)

dom = 20

data = np.loadtxt('/home/nricke/Pictures/tafel/NGCC_O2vsRate_data.csv',delimiter=',') #data[:,0] = log(pO2), data[:,1] = log(J)

J10_data = 10.**data[:,1]
pO2_data = 10.**data[:,0]

print(data)
print(J10_data)
print(pO2_data)

#plt.plot(data[:,0], data[:,1])
#plt.show()


## Linear fit of Log-Log data, as done in the N-GCC paper
"""
regr = linear_model.LinearRegression()
regr.fit(data[:,0].reshape(-1,1), data[:,1])
#regr.fit(data)
print(data[:,0].reshape(-1,1))
print('Coefficients: \n', regr.coef_)
"""

## Fit pO2 vs Rate data for 2-step process without disorder

sig = 0.01
V = -0.5
k_list = [0.001,.01,0.01,0.001] 
dis = IH.DisorderGCC(k_list, sig)
dis.pH = 1

fit = DataFit(dis, V, pO2dom=pO2_data, Jdata=J10_data)

bnds = ((0, None),(0,None),(0,None),(0,None))
res = opt.minimize(fit.evalfit, x0=k_list, bounds=bnds, method='SLSQP')
#res = opt.minimize(fit.evalfit, x0=k_list, bounds=bnds, method='TNC')

print(res.x)
print(res.success)
dis.setConsts(res.x)
print(fit.evalfit(res.x))

O2_dense = np.linspace(0.00001,1.,100)
rate_fO2 = fit.O2vsRate(res.x, pO2=O2_dense)
plt.plot(O2_dense, rate_fO2)


## Linear fit of O2 vs Rate data (on non-log axis). The simplest mechanism would predict a linear dependence
slope, intercept, r_value, p_value, std_err = stats.linregress(pO2_data, J10_data)
line = slope*pO2_data+intercept
plt.plot(pO2_data, line, c='red')



## Fit pO2 vs Rate for same system as above, except with intrinsic disorder in intermediate energies

fitdis = DataFit(dis, V, pO2dom=pO2_data, Jdata=J10_data, disorder=True)
fitdis.emech.sig = 0.003

bnds = ((0, None),(0,None),(0,None),(0,None))
resdis = opt.minimize(fitdis.evalfit, x0=res.x, bounds=bnds, method='SLSQP')

dis.setConsts(resdis.x)
print(resdis.x)
print(resdis.success)
print(fitdis.evalfit(resdis.x))

rate_fO2 = fitdis.O2vsRate(resdis.x, pO2=O2_dense)
plt.plot(O2_dense, rate_fO2, color='green')

rate_fO2 = fitdis.O2vsRate(res.x, pO2=O2_dense)
plt.plot(O2_dense, rate_fO2, color='blue')

## Plot actual data and set plot parameters
plt.scatter(pO2_data, J10_data, c='orange')
#plt.xlim([pO2_data[0]-0.01, pO2_data[-1]+0.01])
plt.show()
