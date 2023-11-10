import Program.Constants as c
from Program.Main import *

rng = np.random.default_rng()

def generateMCArrays(mean, sd):
    return rng.lognormal(np.log(mean), sd/mean, 100)


def MonteCarloErr(f_gas, df_gas, M_200, dM_200, M_hse, dM_hse):
    # Parameters to be Monte-Carlo'd Overwritten

    Omega_b = generateMCArrays(c.Omega_b, c.dOmega_b)
    f_star = generateMCArrays(c.f_star, c.df_star)

    M_hse = generateMCArrays(M_hse, dM_hse)
    f_gas = generateMCArrays(f_gas, df_gas)
    M_200 = generateMCArrays(M_200, dM_200)

    # Now we have overwritten the values to be arrays, we can run the calculation again...

    m_Gas = f_gas * M_hse
    f_ourGas = m_Gas / M_200
    f_b = f_star + f_ourGas
    Omega_m = Omega_b / f_b

    return np.std(Omega_m)


