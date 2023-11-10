import numpy as np
import pandas as pd
from Program.Constants import *
from Program.Queries import *
import matplotlib.pyplot as plt
from Program.Equations import *

dir = "Data/Eckert.csv"
df = pd.read_csv(dir)
Omega_MVals = np.zeros(len(df.index))

ra = ['12 57 11.6', '13 48 50.48', '15 10 56.2', '15 58 14.38', '17 12 15.04',
      '19 20 45.3', '03 42 39.6', '04 31 11.9', '08 17 24.5', '00 41 36.21']
dec = ['-17 24 34', '+26 35 07.4', '+05 44 42', '+27 12 57.8', '+64 03 10.6',
       '+43 57 43', '-53 37 50', '-61 24 23', '-07 30 46', '-09 19 30.4']
size = ['0d23m0s', '0d23m0s', '0d6m0s', '0d16m0s', '0d10m0s', '0d16m0s',
        '0d22m0s', '0d26m0s', '0d10m0s', '0d22m0s']

for i, row in df.iterrows():
    customSimbad = Simbad()
    customSimbad.add_votable_fields('ra(s)', 'dec(s)',
                                    'distance_result', 'otype', 'rv_value', 'z_value')
    customSimbad.remove_votable_fields('coordinates')

    galaxiesInCluster = customSimbad.query_region(coord.SkyCoord(ra[i], dec[i], unit=(u.hourangle, u.deg)),
                                                  radius=size[i]).to_pandas().dropna()

    galaxiesInCluster = galaxiesInCluster[galaxiesInCluster['OTYPE'].isin(acceptableTypes)]
    galaxiesInCluster = galaxiesInCluster[np.abs(galaxiesInCluster['Z_VALUE'] - df['Redshift'][i]) <= 0.035]

    


    M_gas = float(row['f_gas']) * (float(row['Mhse'])) * 1E14 * M_sun
    f_gas = M_gas / M_200(np.std(galaxiesInCluster['RV_VALUE']) * 1000, row['Redshift'])
    f_b = f_star + f_gas
    Omega_m = Omega_b / f_b
    Omega_MVals[i] = Omega_m

    print("\nDATA FOR : " + row['Name'])
    print("--- f_gas: " + str(f_gas))
    print("--- M_tot value: " + str(row['M_tot * 10^14 * M_sun'] * 1E14 * M_sun))
    print("--- omega m : " + str(Omega_m))
    print("-------------------------\n")

    numHistBins = 40

    # Now we generate some graphs of the cluster

    fig2, ax2 = plt.subplots(1, 1)
    ax2.hist(galaxiesInCluster['Z_VALUE'], bins=numHistBins)
    ax2.set_ylabel('num')
    ax2.set_xlabel('Z_VALUE')
    ax2.set_title('Z Values for ' + row['Name'])
    plt.show()

    time.sleep(1.5)

FinalValuesDf = pd.DataFrame(df['Name'].values, Omega_MVals)
print("Mean Omega_M : " + str(np.mean(Omega_MVals)))




