import numpy as np
import pandas as pd
import tafel.DisorderMech as DM


## the image for this data is in the same folder as the read_csv
## x: log(J), y: E/V vs NHE
data = pd.read_csv('~/pymod/tafel/data/Jackson_HER_data.csv', names=['logJ','VvsNHE'])
## this data actually has a 3rd column, which is pH. Need to add that data to this dataframe
## every pH unit has 6 data points, except pH = 3
## the data is sorted in order of ascending pH
## the data also skips pH = 5. Overall there are 65 data points

pH_values = [1,2,3,4,6,7,8,9,10,11,12] #1-12 except 5
pH_list = []
for i in pH_values:
    for j in range(6):
        pH_list.append(i)

pH_list.remove(3)
print(pH_list)
print(len(pH_list))

data = data.assign(pH=pH_list)
print(data[:29])
print(data[29:])
