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

font = {'size':18}
mpl.rc('font',**font)
rcParams['figure.figsize'] = 10,10


label_string = "t,A,B,OO,OA+AO,OB+BO,AA,AB+BA,BB,OOO,OOA+AOO,OOB+BOO,AOA,AOB+BOA,BOB,OAO,OAA+AAO,OAB+BAO,AAA,AAB+BAA,BAB,OBO,OBA+ABO,OBB+BBO,ABA,ABB+BBA,BBB"
labels = label_string.split(',')
print(labels)

indir = "/home/nricke/work/tafel/KMC/500x_rep/"

df_list = []
for infile in os.listdir(indir):
    if infile[:3] == "run" and infile[-3:] == "csv":
        #print(infile)
        df_list.append(pd.read_csv(indir+infile, header=None, names=labels))
        
print(len(df_list))

df = pd.concat(df_list)
#print(df)

df_x = df[["A","B","OO","OA+AO","OB+BO","AA","AB+BA","BB"]].values


df_OAB = df["OAB+BAO"]
df_AAB = df["AAB+BAA"]
df_BAB = df["BAB"]
df_OBA = df["OBA+ABO"]
df_ABA = df["ABA"]
df_ABB = df["ABB+BBA"]
df_BBB = df["BBB"]
df_AAA = df["AAA"]

#poly = skprep.PolynomialFeatures(2)
#X_poly = poly.fit_transform(df_x)
X_poly = df_x
print("type of df_x: ", type(df_x))

#sc_x = skprep
#X = skprep.scale(X_poly)

#y = np.array(df_AAA)
y = np.array(df_OAB)

sc_x = skprep.StandardScaler()
sc_y = skprep.StandardScaler()
X_std = sc_x.fit_transform(X_poly)
#print(X_std)
#y_std = sc_y.fit_transform(y)

print("y: ", y)


#regr = sklm.ElasticNet(alpha=0.05, l1_ratio=0.5, random_state=0)
regr = sklm.LinearRegression()
regr.fit(X_std, y)

y_pred = regr.predict(X_std)

# The coefficients
print('Coefficients: \n', regr.coef_)
# The mean squared error
print("Mean squared error: %.8f" % skm.mean_squared_error(y, y_pred))
# Explained variance score: 1 is perfect prediction
print('Variance score: %.4f' % skm.r2_score(y, y_pred))

plt.scatter(y_pred, y)
#plt.show()

#kf = skms.KFold(n_splits=5)
#results = model_selection.cross_val_score(regr, X, y, cv=kf, scoring="neg_mean_squared_error")
#print("MSE: %.3f (%.3f)") % (results.mean(), results.std())
#
#for train_index, test_index in kf.split(X):
#    X_train, X_test = X[train_index], X[test_index]
#    y_train, y_test = y[train_index], y[test_index]
#    regr = sklm.ElasticNet(alpha=0.05, l1_ratio=0.5, random_state=0)
#    regr.fit(X_train, y_train)
#    y_predict = regr.predict(X_test)

plt.ylabel("OAB from KMC")
plt.xlabel("Fit")
plt.savefig("OAB_o1.png", transparent=True, bbox_inches='tight', pad_inches=0.05)


