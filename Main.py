
from matplotlib import pyplot as plt
from Equations import *
from astroquery.simbad import Simbad
import pandas as pd
import astropy.coordinates as coord
import time
import astropy.units as u

previousQueryTime = time.time()


def findClusters(saveToFile=False, name='clusterQuery.csv'):

    # Setting up query fields required

    customSimbad = Simbad()
    customSimbad.add_votable_fields('z_value', 'ra', 'dec', 'dim_majaxis', 'dim_minaxis')
    customSimbad.remove_votable_fields('coordinates')
    qry = "children > 1000 & otypes in ('ClG')"
    clusters = customSimbad.query_criteria(qry).to_pandas()
    clusters = clusters.dropna()

    if saveToFile:
        clusters.to_csv(name)

    return clusters


def getGalaxiesFromCluster(ra, dec, majaxis, clusterZ):
    global previousQueryTime

    # set up custom query and define which columns we want
    customSimbad = Simbad()
    customSimbad.add_votable_fields('ra(s)', 'dec(s)',
                                    'distance_result', 'otype', 'rv_value', 'z_value')
    customSimbad.remove_votable_fields('coordinates')

    my_radius = '0d' + str(int(majaxis)) + 'm0s'
    print("radius : " + str(my_radius))

    currentQueryTime = time.time()

    if currentQueryTime - previousQueryTime < 1.5:
        print("waiting")
        time.sleep(1.5-currentQueryTime+previousQueryTime)

    df = customSimbad.query_region(coord.SkyCoord(ra, dec, unit=(u.hourangle, u.deg)),
                                             radius=my_radius).to_pandas()
    previousQueryTime = currentQueryTime

    # Clean up data

    # This worked well for the a2029 galaxy
    # return gal_df[(gal_df['Z_VALUE'] < 0.09) * (gal_df['Z_VALUE'] > 0.065)]

    # Were going to cull the data, from our investigations into a2029, we found
    # z +- 0.013 is a good separation for the orbiting galaxies.

    redshiftDifferenceCutoff = 0.013

    print(pd.Series({c: df[c].unique() for c in df})['OTYPE'])

    gal_df = df[df['OTYPE'] == 'Galaxy']
    gal_df = gal_df.dropna()
    gal_df = gal_df[np.abs(gal_df['Z_VALUE'] - clusterZ) <= redshiftDifferenceCutoff]

    return gal_df


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


c = findClusters(saveToFile=True)

clustersSet = pd.read_csv('clusterQuery.csv')

for i in clustersSet.index:
    clusterData = getGalaxiesFromCluster(clustersSet['RA'][i],
                                         clustersSet['DEC'][i],
                                         clustersSet['GALDIM_MAJAXIS'][i],
                                         clustersSet['Z_VALUE'][i])

    # if len(clusterData.index) < 10:
    #     print("\nInsufficient Number of galaxies in cluster : " + str(clustersSet['MAIN_ID'][i]))
    #     continue


    print(clusterData)
    CalculateOmegaM(clusterData, nameOfGalaxy=clustersSet['MAIN_ID'][i], makeGraphs=True)
    break


print("\nFinished Queries")
