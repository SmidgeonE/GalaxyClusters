# %%

import Equations
from Main import *

rng = np.random.default_rng()


def generateMCArrays(mean, sd):
    return rng.lognormal(np.log(mean), sd/mean, 100)


# Parameters to be Monte-Carlo'd Overwritten

Equations.G = generateMCArrays(Equations.G, Equations.dG)
Equations.H_0 = generateMCArrays(Equations.H_0, Equations.dH_0)
Equations.Omega_M = generateMCArrays(Equations.Omega_M, Equations.dOmega_M)
Equations.Omega_L = 1-Omega_M
Equations.h = H_0 / 100
Equations.Omega_b = generateMCArrays(Equations.Omega_b, Equations.dOmega_b)
Equations.f_gas = generateMCArrays(Equations.f_gas, Equations.df_gas)
Equations.f_star = generateMCArrays(Equations.f_star, Equations.df_star)
Equations.M_hse = generateMCArrays(Equations.M_hse, Equations.dM_hse)
Equations.M_gas = f_gas * M_hse

# Now we have overwritten the values to be arrays, we can run the calculation again..

mVals = CalculateOmegaM()

print("Monte Carlo Mean for Omega M : " + str(np.mean(mVals)))
print("Monte Carlo Error for Omega M : " + str(np.std(mVals)))
