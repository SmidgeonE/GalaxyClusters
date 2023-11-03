import pandas as pd
from scipy import optimize
from Program.Equations import *


dataDir = "Data/Mantzetal.csv"
data = pd.read_csv(dataDir)

# Removing the asterisked values, they are described in the paper as being inaccurate.

for i, row in data.iterrows():
    if row['Name'][-1] == '∗':
        data = data[data.index != i]


def f_gasFunc(zM, f_gas0, s, k):
    z, M = zM
    return f_gas0 * np.power(1 + z, -s) * np.power(M/(10E15 * (1/h) * M_sun), k)


x = optimize.curve_fit(f_gasFunc, [data['Redshift'], np.array(data['M_tot * 10^14 * M_sun']) * 10E14 * M_sun],
                   data['f_gas'])[0]


# These constants seem to be correct, the f_gas0 is 0.14 where Majmudar uses 0.15.

f_gas0 = x[0]
s = x[1]
k = x[2]

