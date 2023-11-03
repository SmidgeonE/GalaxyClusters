import numpy as np
import Program.fGas as fGas

G = 6.67408E-11
H_0 = 67.80
Omega_M = 0.15
Omega_L = 1-Omega_M
h = H_0 / 100
M_sun = 2E30

# Parameters to be investigated for procedural generation
# Although omega_b shouldnt vary much, implied by lambdaCDM

Omega_b = 0.0125 / h**2
f_gas = 0.152
f_star = 0.015


# parameters from paper

f_paperGas = 0.152
M_hse = 12.25E14 * 2E30

# Errors

dG = 0.0005 * G
dH_0 = 0.77
dOmega_M = 0.04
dOmega_L = dOmega_M
dOmega_b = 0.001
df_gas = 0.006
df_star = 0.005
dM_hse = 0.49E14 * 2E30


def convertDist(x, unit, to):
    valInM = 0

    if unit == "km":
        valInM = x * 1000
    elif unit == "Mpc":
        valInM = 3.08568E22 * x

    # Returns

    if to == "km":
        return valInM * 1000
    elif to == "Mpc":
        return valInM * 3.08568E22 * valInM


def F_gas(z, M):
    return fGas.f_gasFunc((z, M), fGas.f_gas0, fGas.s, fGas.k)


def gauss(data, maxCount, range):
    data = data.to_numpy()
    mean = np.mean(data)
    sd = np.std(data)

    return maxCount * np.exp(-0.5 * (range-mean)**2 / sd**2), sd


def M_200(sigma_v, z):
    return np.sqrt(81/(800*G**3*np.pi*rho_z(z))) * sigma_v**3


def R_200(sigma_v, z):
    return np.cbrt(3*M_200(sigma_v, z) / (800 * np.pi * rho_z(z)))


def rho_z(z):
    # print("rho_C : " + str(3*(H(z)**2) / (8*np.pi*G)))
    return 3*(H(z)**2) / (8*np.pi*G)


def H(z):
    # print("H_z :" + str(H_0 * np.sqrt(Omega_M * (1+z)**3 + Omega_L) * 3.08568E22 * 1000))
    return H_0 * np.sqrt(Omega_M * (1+z)**3 + Omega_L) / 3.08568E22 * 1000


