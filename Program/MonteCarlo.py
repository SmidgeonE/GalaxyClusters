from Program.Constants import f_gas, M_hse, dG, dH_0, dOmega_M, dOmega_b, df_gas, \
    df_star, dM_hse
from Main import *

rng = np.random.default_rng()


def generateMCArrays(mean, sd):
    return rng.lognormal(np.log(mean), sd/mean, 100)


# Parameters to be Monte-Carlo'd Overwritten

G = generateMCArrays(G, dG)
H_0 = generateMCArrays(H_0, dH_0)
Omega_M = generateMCArrays(Omega_M, dOmega_M)
Omega_L = 1 - Omega_M
h = H_0 / 100
Omega_b = generateMCArrays(Omega_b, dOmega_b)
f_gas = generateMCArrays(f_gas, df_gas)
f_star = generateMCArrays(f_star, df_star)
M_hse = generateMCArrays(M_hse, dM_hse)
M_gas = f_gas * M_hse

# Now we have overwritten the values to be arrays, we can run the calculation again..

mVals = CalculateAndPlotOmegaM()

print("Monte Carlo Mean for Omega M : " + str(np.mean(mVals)))
print("Monte Carlo Error for Omega M : " + str(np.std(mVals)))
