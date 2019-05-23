import pandas as pd
import numpy as np

#import keras
from keras.models import Sequential
from keras.layers import Dense
from keras.wrappers.scikit_learn import KerasRegressor

from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

rn_seed = 74
np.random.seed(rn_seed)

#infile = sys.argv[1]
#infile = "/home/nricke/work/tafel/KMC/k_run.csv"
infile = "/home/nricke/work/tafel/KMC/k_run_prep.csv"
df = pd.read_csv(infile, index_col="Unnamed: 0")

X = df[["OA+AO","OB+BO","AA","AB+BA","BB", "B3_Estimate"]].values
Y = df["OAB+BAO"].values

def baseline_model():
    ffnn = Sequential()
    ffnn.add(Dense(36, input_dim=6, activation='relu'))
    ffnn.add(Dense(36, input_dim=36, activation='relu'))
    ffnn.add(Dense(18, input_dim=36, activation='relu'))
    ffnn.add(Dense(1, input_dim=18, activation='linear'))
    ffnn.compile(loss='mean_squared_error', optimizer='adam', metrics=[])
    return ffnn

"""
# train test split
X_train, X_test, y_train, y_test = skms.train_test_split(X, Y, test_size=0.2, random_state=1)
ffnn.fit(X_train, Y_train, epochs=1000, batch_size=500)
scores = model.evaluate(X_train, Y_train)
predictions = model.predict(X_test)
"""

estimators = []
estimators.append(('standardize', StandardScaler()))
estimators.append(('mlp', KerasRegressor(build_fn=baseline_model, epochs=5, batch_size=50, verbose=0)))
pipeline = Pipeline(estimators)
kfold = KFold(n_splits=3, random_state=rn_seed)
results = cross_val_score(pipeline, X, Y, cv=kfold)
print("Results: %.2f (%.2f) MSE" % (results.mean(), results.std()))




print()
