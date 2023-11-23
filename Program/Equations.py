import numpy as np
import Program.fGas as fGas
from Program.Constants import G, H_0, Omega_M, Omega_L


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


def M_200(sigma_v, z, omega_M=Omega_M):
    return np.sqrt(81 / (800 * G ** 3 * np.pi * rho_z(z, omega_M))) * sigma_v**3


def R_200(sigma_v, z, omega_M=Omega_M):
    return np.cbrt(3*M_200(sigma_v, z, omega_M) / (800 * np.pi * rho_z(z, omega_M)))


def rho_z(z, omega_M):
    # print("rho_C : " + str(3*(H(z)**2) / (8*np.pi*G)))
    return 3*(H(z, omega_M)**2) / (8 * np.pi * G)


def H(z, omega_M):
    return H_0 * np.sqrt(omega_M * (1 + z) ** 3 + (1-omega_M)) / 3.08568E22 * 1000


