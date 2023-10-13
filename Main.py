# %%

from matplotlib import pyplot as plt
from Equations import *
from astroquery.simbad import Simbad
import astropy.coordinates as coord


def findClusters(saveToFile=False, name='clusterQuery.csv'):

    # Querying to find galaxy clusters

    customSimbad = Simbad()

    qry = ("region(circle, 29.20 -0.214, 0.15d) &"
           " otypes in ('ClG', 'C?G')")
    clusters = customSimbad.query_criteria(qry).to_pandas()
    clusters = clusters.loc[:, ['MAIN_ID', 'RA', 'DEC']]

    print(clusters)

    if saveToFile:
        clusters.to_csv(name)

    return clusters


def getClusterData(ra=227.745, dec=5.761):
    # set up custom query and define which columns we want
    customSimbad = Simbad()
    customSimbad.add_votable_fields('ra(d)', 'dec(d)', 'distance_result', 'otype', 'rv_value', 'z_value')
    customSimbad.remove_votable_fields('coordinates')

    # use 20 arcmin as an example
    my_radius = '0d20m0s'

    # run query
    result_table = customSimbad.query_region(coord.SkyCoord(ra, dec, unit='deg'), radius=my_radius)
    df = result_table.to_pandas()

    # Select only galaxies
    gal_df = df[df['OTYPE'] == 'Galaxy']
    return gal_df[(gal_df['Z_VALUE'] < 0.09) * (gal_df['Z_VALUE'] > 0.065)]


cluster = getClusterData()


def CalculateOmegaM(makeGraphs=False):

    # Now we generate some graphs of the cluster

    fig2, ax2 = plt.subplots(1, 1)
    ax2.hist(cluster['Z_VALUE'], bins=100)
    ax2.set_ylabel('num')
    ax2.set_xlabel('Z_VALUE')
    if makeGraphs: plt.show()

    # Now we plot a gaussian on the histogram for recession vel

    recessionVel = cluster['RV_VALUE']
    fig, ax = plt.subplots(1, 1)
    histData = ax.hist(recessionVel, bins=100)
    ax.set_ylabel('num')
    ax.set_xlabel('RV_VALUE')

    velocityRange = np.linspace(np.min(recessionVel), np.max(recessionVel), len(recessionVel))
    gaussVals, sd = gauss(recessionVel, np.max(histData[0]), velocityRange)
    ax.plot(velocityRange, gaussVals)
    if makeGraphs: plt.show()

    averageZ = cluster['Z_VALUE'].mean()
    M_200val = M_200(sd * 1000, averageZ)
    f_ourGas = M_gas / M_200val
    f_b = f_star + f_ourGas
    Omega_m = Omega_b / f_b

    print("M_200 value: " + str(M_200val))
    print("omega m : " + str(Omega_m))

    return Omega_m


# %%

c = findClusters(saveToFile=True)
