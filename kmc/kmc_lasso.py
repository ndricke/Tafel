import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from pylab import rcParams

import sklearn.preprocessing as skprep
import sklearn.linear_model as sklm
import sklearn.model_selection as skms
import sklearn.metrics as skm
from sklearn.decomposition import PCA

font = {'size':18}
mpl.rc('font',**font)
rcParams['figure.figsize'] = 10,10

"""
label_string = "t,A,B,OO,OA+AO,OB+BO,AA,AB+BA,BB,OOO,OOA+AOO,OOB+BOO,AOA,AOB+BOA,BOB,OAO,OAA+AAO,OAB+BAO,AAA,AAB+BAA,BAB,OBO,OBA+ABO,OBB+BBO,ABA,ABB+BBA,BBB"
labels = label_string.split(',')
print(labels)

#indir = "/home/nricke/work/tafel/KMC/500x_rep/"
indir = "/home/nricke/work/tafel/KMC/LV/k_run"

df_list = []
#for infile in os.listdir(indir):
for path, subdirs, files in os.walk(indir):
    print(path, subdirs, files)
    for name in files:
        if name[:3] == "run" and name[-3:] == "csv":
            df_list.append(pd.read_csv(path+'/'+name, header=None, names=labels))

print(len(df_list))

df = pd.concat(df_list)
df.to_csv("k_run.csv")
sys.exit(-1)
"""

ep = 10.**-10

df = pd.read_csv("k_run.csv")
B3_string = sys.argv[1] # the 3-site correlation term to estimate



#df_x = df[["A","B","OO","OA+AO","OB+BO","AA","AB+BA","BB"]].values

#df_x = df[["A","B","OA+AO","OB+BO","AA","AB+BA","BB"]]
#df_x["A"] = 1./(df_x["A"]+ep)
#df_x["B"] = 1./(df_x["B"]+ep)
#df_x["O"] = 1./(df_x["O"]+ep)



df["B3_Estimate"] = df[B3_left2]*df[B3_right2]/df[B3_middle]

df_x = df[["OA+AO","OB+BO","AA","AB+BA","BB"]].values

#df = df[["OA+AO","OB+BO","AA","AB+BA","BB","OAB+BAO","AAB+BAA","BAB","OBA+ABO","ABA","ABB+BBA","BBB","AAA"]]

#df_AAA = df["OAB+BAO"]
#df_AAA = df["AAB+BAA"]
#df_AAA = df["BAB"]
#df_AAA = df["OBA+ABO"]
#df_AAA = df["ABA"]
#df_AAA = df["ABB+BBA"]
#df_AAA = df["BBB"]
#df_AAA = df["AAA"]
df_B3 = df[B3]


poly = skprep.PolynomialFeatures(3)
X_poly = poly.fit_transform(df_x)
#X_poly = df_x
print("type of df_x: ", type(df_x))

#sc_x = skprep
#X = skprep.scale(X_poly)

#y = np.array(df_AAA)
y = np.array(df_AAA)

sc_x = skprep.StandardScaler()
sc_y = skprep.StandardScaler()
X_std = sc_x.fit_transform(X_poly)
#print(X_std)
#y_std = sc_y.fit_transform(y)

print("y: ", y)

#pca = PCA(0.9999).fit(X_std)
#plt.plot(np.cumsum(pca.explained_variance_ratio_))
#plt.xlabel("number of components")
#plt.ylabel("cumulative explained variance")
#plt.show()

#X_std = pca.transform(X_std)
#print(X_std.shape)

#alpha_list = [10.**-4, 5*10.**-5, 10.**-5]
#for al in alpha_list:
#regr = sklm.ElasticNet(alpha=al, l1_ratio=0.9, random_state=0)
regr = sklm.Lasso(alpha=al, random_state=0, max_iter=2000)
#regr = sklm.LinearRegression()

regr.fit(X_std, y)
y_pred = regr.predict(X_std)
# The coefficients
print('Coefficients: \n', regr.coef_)
# The mean squared error
print("Mean squared error: %.8f" % skm.mean_squared_error(y, y_pred))
# Explained variance score: 1 is perfect prediction
print('Variance score: %.4f' % skm.r2_score(y, y_pred))

#kf = skms.KFold(n_splits=5)
#results = skms.cross_val_score(regr, X_std, y, cv=kf, scoring="neg_mean_squared_error")
#print(al, np.mean(results), np.std(results))

#plt.scatter(y_pred, y)
#plt.ylabel("OAB from KMC")
#plt.xlabel("Fit")
#plt.xlim([-0.004,np.max(y_pred)+0.005])
#plt.ylim([-0.004, np.max(y)+0.005])
#plt.show()
#plt.savefig("AAB_o3lasso_krun.png", transparent=True, bbox_inches='tight', pad_inches=0.05)



#"""
