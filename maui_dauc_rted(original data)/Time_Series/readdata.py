import pandas as pd

df = pd.read_csv (r'TIMESERIES_RT_LOAD.csv')

# import csv

# file = open('TIMESERIES_RT_LOAD.csv')

# csvreader = csv.reader(file)

# header = []
# header = next(csvreader)
# for row in csvreader:
#     print(row)

# get bus name
df.columns=df.columns[2:]