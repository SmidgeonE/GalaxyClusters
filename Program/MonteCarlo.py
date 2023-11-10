from Constants import *
from Main import *

rng = np.random.default_rng()

def generateMCArrays(mean, sd):
    return rng.lognormal(np.log(mean), sd/mean, 100)


def MonteCarloErr(self, f_gas, df_gas, M_200, dM_200, M_hse, dM_hse):
    # Parameters to be Monte-Carlo'd Overwritten

    self.Omega_b = generateMCArrays(self.Omega_b, dOmega_b)
    self.f_star = generateMCArrays(self.f_star, df_star)

    M_hse = generateMCArrays(M_hse, dM_hse)
    f_gas = generateMCArrays(f_gas, df_gas)
    M_200 = generateMCArrays(M_200, dM_200)

    # Now we have overwritten the values to be arrays, we can run the calculation again...

    m_Gas = f_gas * M_hse
    f_ourGas = m_Gas / M_200
    f_b = f_star + f_ourGas
    Omega_m = Omega_b / f_b


    print("Monte Carlo Mean for Omega M : " + str(np.mean(Omega_m)))
    print("Monte Carlo Error for Omega M : " + str(np.std(Omega_m)))
