import numpy as np

G = 6.67408E-11
H_0 = 100
Omega_M = -0.15
Omega_L = -Omega_M

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

def M(sigma_v, R):
    return (3/G) * sigma_v**2 * R

def R():
    return 0


def rho_z(z):
    return 3*H(z)**2 / 8*np.pi*G

def H(z):
    return H_0 * np.sqrt(Omega_M * (1+z)**3 + Omega_L)


