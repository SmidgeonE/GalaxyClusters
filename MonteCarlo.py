# %%

import numpy as np
import Equations
from Main import *

rng = np.random.default_rng()


def generateMCArrays(mean, sd):
    return rng.lognormal(np.log(mean), sd/mean, 100)


# Parameters to be Monte-Carlo'd

G = generateMCArrays(Equations.G, Equations.dG)
H_0 = generateMCArrays(Equations.H_0, Equations.dH_0)
Omega_M = generateMCArrays(Equations.Omega_M, Equations.dOmega_M)
Omega_L = 1-Omega_M
h = H_0 / 100
Omega_b = generateMCArrays(Equations.Omega_b, Equations.dOmega_b)
f_gas = generateMCArrays(Equations.f_gas, Equations.df_gas)
f_star = generateMCArrays(Equations.f_star, Equations.df_star)
M_hse = generateMCArrays(Equations.M_hse, Equations.dM_hse)
M_gas = f_gas * M_hse

# Now we have overwritten the values to be arrays, we can run the calculation again..

getDataAndCalculateOmegaM()
