# %%

import numpy as np
from matplotlib import pyplot as plt
from Equations import *
from astropy.cosmology import WMAP9 as cosmo
from astroquery.simbad import Simbad
import astropy.coordinates as coord

# set up custom query and define which columns we want
customSimbad = Simbad()
customSimbad.add_votable_fields('ra(d)', 'dec(d)', 'distance_result', 'otype', 'rv_value', 'z_value')
customSimbad.remove_votable_fields('coordinates')

# set RA, Dec and search radius
my_ra = 227.745
my_dec = 5.761

# use 20 arcmin as an example
my_radius = '0d20m0s'

# run query
result_table = customSimbad.query_region(coord.SkyCoord(my_ra, my_dec, unit='deg'), radius=my_radius)
df = result_table.to_pandas()

# print(df)
print("The columns in the data frame are:\n",df.columns.values,"\n")

# Select only galaxies
gal_df = df[df['OTYPE'] == 'Galaxy']
print('Selected',len(gal_df),'galaxies\n')

cluster = gal_df[(gal_df['Z_VALUE'] < 0.09) * (gal_df['Z_VALUE'] > 0.065)]

plt.show()

fig2, ax2 = plt.subplots(1,1)
ax2.hist(cluster['Z_VALUE'], bins=100)
ax2.set_ylabel('num')
ax2.set_xlabel('Z_VALUE')
plt.show()

# Now we plot a gaussian on the histrogram for recession vel

recessionVel = cluster['RV_VALUE']
fig, ax = plt.subplots(1,1)
histData = ax.hist(recessionVel, bins=100)
ax.set_ylabel('num')
ax.set_xlabel('RV_VALUE')

range = np.linspace(np.min(recessionVel), np.max(recessionVel), len(recessionVel))
gaussVals, sd = gauss(recessionVel, np.max(histData[0]), range)
ax.plot(range, gaussVals)
plt.show()


averageZ = cluster['Z_VALUE'].mean()

M_200 = M_200(sd * 1000, averageZ)
print("M_200 value: " + str(M_200))

f_ourGas = M_gas / M_200
f_star = 0.015
f_b = f_star + f_ourGas


print(f_ourGas)
print("omega m : " + str(Omega_b / f_b))


