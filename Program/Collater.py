# Collating the query objects for clusters
# Because large querys dont return all valid objects, for some reason
# Ok I thought some data sets returned values that others didn't, turns out the smaller sets are just subsets
# of the largest set.

import pandas as pd


dirs = ["Data/clusterQuery20children.csv",
        "Data/clusterQuery21children.csv",
        "Data/clusterQuery40children.csv",
        "Data/clusterQuery50children.csv",
        "Data/clusterQuery100children.csv"]


df = pd.read_csv(dirs[0])

for d in dirs:
    if d == dirs[0]:
        continue

    df = pd.concat([df, pd.read_csv(d)])

df = df.drop_duplicates(subset=["MAIN_ID"])
print(len(df[df['MAIN_ID'] == "NAME Sausage Cluster"].index))
print(len(df.index))



