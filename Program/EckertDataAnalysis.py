import time

import numpy as np
import pandas as pd
from Program.Constants import *
from Program.Queries import *
import matplotlib.pyplot as plt
from Program.Equations import *
from Program.bootstrap import *
from Program.MonteCarlo import *
import scienceplots
import os

dir = "Data/Eckert.csv"
df = pd.read_csv(dir)
Omega_MVals = np.zeros(len(df.index))
MCErrors = np.zeros(len(df.index))

ra = ['12 57 11.6', '13 48 50.48', '15 10 56.2', '15 58 14.38', '17 12 15.04',
      '19 20 45.3', '03 42 39.6', '04 31 11.9', '08 17 24.5', '00 41 36.21']
dec = ['-17 24 34', '+26 35 07.4', '+05 44 42', '+27 12 57.8', '+64 03 10.6',
       '+43 57 43', '-53 37 50', '-61 24 23', '-07 30 46', '-09 19 30.4']
size = ['0d23m0s', '0d23m0s', '0d6m0s', '0d16m0s', '0d10m0s', '0d16m0s',
        '0d22m0s', '0d26m0s', '0d30m0s', '0d22m0s']
cutoffUpper = [0.055, 0.0675, 0.0875, 0.1, 0.0875, 0.065, 0.0675, 0.07, 0.009, 0.0625]
cutoffLower = [0.0325, 0.055, 0.07, 0.08, 0.0725, 0.0425, 0.05, 0.05, 0.05, 0.045]


# This is the very first Eckert data analysis, this gave very low values.

# f_gas = df['f_gas'].values * df['Mhse'].values / df['M_tot * 10^14 * M_sun'].values
# f_b = f_star + f_gas
# Omega_m = Omega_b / f_b
# print(np.std(Omega_m))


def AnalyseData():
    print(Omega_M)
    for i, row in df.iterrows():
        clusterDir = "Data/GalaxySets/" + str(row['Name']) + ".csv"

        if not os.path.exists(clusterDir):
            customSimbad = Simbad()
            customSimbad.add_votable_fields('ra(s)', 'dec(s)',
                                            'distance_result', 'otype', 'rv_value', 'z_value')
            customSimbad.remove_votable_fields('coordinates')

            try:
                galaxiesInCluster = customSimbad.query_region(coord.SkyCoord(ra[i], dec[i], unit=(u.hourangle, u.deg)),
                                                              radius=size[i]).to_pandas().dropna()
            except:
                print("     ### HTTPS too many requests. Waiting")
                time.sleep(20)
                galaxiesInCluster = customSimbad.query_region(coord.SkyCoord(ra[i], dec[i], unit=(u.hourangle, u.deg)),
                                                              radius=size[i]).to_pandas().dropna()

            galaxiesInCluster.to_csv("Data/GalaxySets/" + str(row['Name']) + ".csv")
            time.sleep(2)
        else:
            print("Found cluster data")
            galaxiesInCluster = pd.read_csv(clusterDir)

        galaxiesInCluster = galaxiesInCluster[galaxiesInCluster['OTYPE'].isin(acceptableTypes)]
        galaxiesInCluster = galaxiesInCluster[(galaxiesInCluster['Z_VALUE'] < cutoffUpper[i])
                                              * (galaxiesInCluster['Z_VALUE'] > cutoffLower[i])]

        M_gas = float(row['f_gas']) * (float(row['Mhse'])) * 1E14 * M_sun


        m_200Val = M_200(np.std(galaxiesInCluster['RV_VALUE']) * 1000, row['Redshift'])
        m_200Err = M_200Err(row['Redshift'], np.std(galaxiesInCluster['RV_VALUE']), BootstrapErr(galaxiesInCluster))

        f_gas = M_gas / m_200Val
        f_b = f_star + f_gas
        Omega_m = Omega_b / f_b
        Omega_MVals[i] = Omega_m
        MCErrors[i] = OmegaMErrorEckert(f_gas, row['f_gasErr'],
                                        m_200Val, m_200Err,
                                        row['Mhse'] * 1E14 * M_sun, row['MhseErr'] * 1E14 * M_sun)

        print("\nDATA FOR : " + row['Name'])
        print("--- f_gas: " + str(f_gas))
        print("--- M_tot value: " + str(row['M_tot * 10^14 * M_sun'] * 1E14 * M_sun))
        print("--- sigma_v error: " + str(BootstrapErr(galaxiesInCluster)))
        print("--- omega m : " + str(Omega_m))
        print("\n--- Monte Carlo For Omega_M: " + str(MCErrors[i]))
        print("-------------------------\n")

        numHistBins = 20

        # Now we generate some graphs of the cluster

        fig2, ax2 = plt.subplots(1, 1)
        ax2.hist(galaxiesInCluster['Z_VALUE'], bins=numHistBins)
        ax2.set_ylabel('num')
        ax2.set_xlabel('Z_VALUE')
        ax2.set_title('Z Values for ' + row['Name'])
        # plt.show()


    meanOmegaM = np.nanmean(Omega_MVals)
    stdOmegaM = np.nanstd(Omega_MVals)
    mcError = np.nanmean(MCErrors)
    print("Mean Omega_M : " + str(meanOmegaM))
    print("Standard Dev of Omega_M :" + str(stdOmegaM))
    print("Monte Carlo Error for Omega_M : " + str(mcError))

    return meanOmegaM, stdOmegaM, mcError


rangeOfOmegaM = np.linspace(0.05, 0.9, 10)
results = np.zeros((10, 3))

for i, val in enumerate(rangeOfOmegaM):
    global Omega_M

    Omega_M = val
    results[i, :] = AnalyseData()


pd.DataFrame(results).to_csv("Data/variedOmegaM.csv")


plt.errorbar(results[:, 0], rangeOfOmegaM, xerr=results[:, 1])




