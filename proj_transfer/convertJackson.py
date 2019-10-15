import pandas as pd

"""
I didn't get the data directly from Megan; she just shared her powerpoint slides, and I took the xy values off them
I needed to add the corresponding pH values using this script
Note that pH=3 only has 5 data points, so I have to remove one
"""

# x: log(J), y: E/V vs NHE
data = pd.read_csv('data_csv/Jackson_HER_data.csv', names=['logJ','VvsNHE'])

pH_values = [1,2,3,4,6,7,8,9,10,11,12] #1-12 except 5
pH_list = []
for i in pH_values:
    for j in range(6):
        pH_list.append(i)
pH_list.remove(3) #pH=3 only has 5 data points

data = data.assign(pH=pH_list)
data.to_csv("Jackson_HERfull_data.csv", index=False)
