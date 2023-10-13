import numpy as np

G = 6.67408E-11
H_0 = 70
Omega_M = 0.15
Omega_L = 1-Omega_M
h = H_0 / 100
Omega_b = 0.0125 / h**2


# TO DO: FIX THESE CONSTS

f_gas = 0.152
f_star = 1
M_hse = 12.25E14 * 2E30
M_gas = f_gas * M_hse




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


def gauss(data, maxCount, range):
    data = data.to_numpy()
    mean = np.mean(data)
    sd = np.std(data)

    return maxCount * np.exp(-0.5 * (range-mean)**2 / sd**2), sd



def M_200(sigma_v, z):
    return np.sqrt(81/(800*G**3*np.pi*rho_z(z))) * sigma_v**3

def rho_z(z):
    # print("rho_C : " + str(3*(H(z)**2) / (8*np.pi*G)))
    return 3*(H(z)**2) / (8*np.pi*G)

def H(z):
    # print("H_z :" + str(H_0 * np.sqrt(Omega_M * (1+z)**3 + Omega_L) * 3.08568E22 * 1000))
    return H_0 * np.sqrt(Omega_M * (1+z)**3 + Omega_L) / 3.08568E22 * 1000


