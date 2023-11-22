import pandas as pd
from scipy import optimize
from Program.Constants import h, M_sun
import matplotlib.pyplot as plt
import numpy as np


# dataDir = "Data/Mantzetal.csv"
dataDir = "Data/Eckert.csv"

data = pd.read_csv(dataDir)

# Removing the asterisked values, they are described in the paper as being inaccurate.

for i, row in data.iterrows():
    if row['Name'][-1] == '∗':
        data = data[data.index != i]


def f_gasFunc(zM, f_gas0, s, k):
    z, M = zM
    return f_gas0 * np.power(1 + z, -s) * np.power(M/(10E15 * (1/h) * M_sun), k)


x = optimize.curve_fit(f_gasFunc,
                       [data['Redshift'], np.array(data['M_tot * 10^14 * M_sun']) * 10E14 * M_sun],
                       data['f_gas'])[0]

# These constants seem to be correct, the f_gas0 is 0.14 where Majmudar uses 0.15.
# Mantz's data is much worse than Eckerts

f_gas0 = x[0]
s = x[1]
k = x[2]


def plotParameters():

    fig, ax = plt.subplots(1, 1)
    massRange = np.linspace(0.1, 7) * 10E45
    zVals = np.linspace(0.03, 0.12, 5)

    for val in zVals:
        ax.plot(massRange, f_gasFunc((val, massRange), f_gas0, s, k))

    ax.set_title("f_gas against mass of cluster from parameterised equation.")
    ax.set_ylabel("f_gas")
    ax.set_xlabel("mass")

    ax.scatter(data['M_tot * 10^14 * M_sun'] * 10E14 * M_sun, data['f_gas'])

    plt.show()



