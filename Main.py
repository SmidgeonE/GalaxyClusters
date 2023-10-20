
from matplotlib import pyplot as plt
from Equations import *
from astroquery.simbad import Simbad
import pandas as pd
import astropy.coordinates as coord
import time
import astropy.units as u

previousQueryTime = time.time()


def findClusters(saveToFile=False, name='clusterQuery.csv'):

    # Querying to find galaxy clusters

    customSimbad = Simbad()
    customSimbad.add_votable_fields('z_value', 'ra', 'dec')
    customSimbad.remove_votable_fields('coordinates')

    qry = ("region(circle, 29.20 -0.214, 0.9d) &"
           " otypes in ('ClG')")
    clusters = customSimbad.query_criteria(qry).to_pandas()
    clusters = clusters.drop(['SCRIPT_NUMBER_ID'], axis=1)
    clusters = clusters[clusters['Z_VALUE'] < 0.15]


    print(clusters)

    if saveToFile:
        clusters.to_csv(name)

    return clusters


def getGalaxiesFromCluster(ra, dec):
    global previousQueryTime

    # set up custom query and define which columns we want
    customSimbad = Simbad()
    customSimbad.add_votable_fields('ra(s)', 'dec(s)',
                                    'distance_result', 'otype', 'rv_value', 'z_value')
    customSimbad.remove_votable_fields('coordinates')

    my_radius = '0d2m0s'

    currentQueryTime = time.time()

    if currentQueryTime - previousQueryTime < 0.5:
        print("waiting")
        time.sleep(0.5-currentQueryTime+previousQueryTime)

    result_table = customSimbad.query_region(coord.SkyCoord(ra, dec, unit=(u.hourangle, u.deg)), radius=my_radius)
    previousQueryTime = currentQueryTime

    df = result_table.to_pandas()

    # Select only galaxies
    gal_df = df[df['OTYPE'] == 'Galaxy']

    return gal_df

    # This worked well for the a2029 galaxy
    # return gal_df[(gal_df['Z_VALUE'] < 0.09) * (gal_df['Z_VALUE'] > 0.065)]


def CalculateOmegaM(galaxiesOfCluster, nameOfGalaxy, makeGraphs=False):

    # Now we generate some graphs of the cluster

    fig2, ax2 = plt.subplots(1, 1)
    ax2.hist(galaxiesOfCluster['Z_VALUE'], bins=100)
    ax2.set_ylabel('num')
    ax2.set_xlabel('Z_VALUE')
    ax2.set_title('Z Values for ' + nameOfGalaxy)
    if makeGraphs:
        plt.show()

    # Now we plot a gaussian on the histogram for recession vel

    recessionVel = galaxiesOfCluster['RV_VALUE']
    fig, ax = plt.subplots(1, 1)
    histData = ax.hist(recessionVel, bins=100)
    ax.set_ylabel('num')
    ax.set_xlabel('RV_VALUE')

    velocityRange = np.linspace(np.min(recessionVel), np.max(recessionVel), len(recessionVel))
    gaussVals, sd = gauss(recessionVel, np.max(histData[0]), velocityRange)
    ax.plot(velocityRange, gaussVals)
    if makeGraphs:
        plt.show()

    averageZ = galaxiesOfCluster['Z_VALUE'].mean()
    M_200val = M_200(sd * 1000, averageZ)
    f_ourGas = M_gas / M_200val
    f_b = f_star + f_ourGas
    Omega_m = Omega_b / f_b

    print("\nDATA FOR : " + nameOfGalaxy)
    print("--- M_200 value: " + str(M_200val))
    print("--- omega m : " + str(Omega_m))

    return Omega_m


# c = findClusters(saveToFile=True)

clustersSet = pd.read_csv('clusterQuery.csv')

for i in clustersSet.index:
    clusterData = getGalaxiesFromCluster(clustersSet['RA'][i], clustersSet['DEC'][i])

    # if len(clusterData.index) < 10:
    #     print("\nInsufficient Number of galaxies in cluster : " + str(clustersSet['MAIN_ID'][i]))
    #     continue


    print(clusterData)
    CalculateOmegaM(clusterData, nameOfGalaxy=clustersSet['MAIN_ID'][i], makeGraphs=True)


print("\nFinished Queries")
