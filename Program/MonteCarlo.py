import numpy as np

import Program.Constants as c
from Program.Equations import *

rng = np.random.default_rng()


def generateMCArrays(mean, sd):
    return rng.lognormal(np.log(mean), sd/mean, 1000)


def OmegaMErrorEckert(f_gas, df_gas, M_200, dM_200, M_hse, dM_hse):
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


def OmegaMErrorGeneral(f_gas, df_gas):
    # Parameters to be Monte-Carlo'd Overwritten

    Omega_b = generateMCArrays(c.Omega_b, c.dOmega_b)
    f_star = generateMCArrays(c.f_star, c.df_star)

    f_gas = generateMCArrays(f_gas, df_gas)

    # Now we have overwritten the values to be arrays, we can run the calculation again...

    f_b = f_star + f_gas
    Omega_m = Omega_b / f_b

    return np.std(Omega_m)


def M_200Err(z, sigma_v, dsigma_v):
    sigma_v = generateMCArrays(sigma_v, dsigma_v)

    mcVals = M_200(sigma_v, z)

    return np.std(mcVals)


def f_gasErr(z, M_200, dM_200):
    M_200 = generateMCArrays(M_200, dM_200)

    mcVals = F_gas(z, M_200)

    return np.std(mcVals)



