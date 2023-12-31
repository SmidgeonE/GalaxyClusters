import pandas as pd
from Program.Equations import *
from astroquery.simbad import Simbad
import astropy.coordinates as coord
import time
import astropy.units as u
import os

previousQueryTime = time.time()
acceptableTypes = ['GinPair', 'Galaxy', 'Compact_Gr_G', 'GroupG', 'LowSurfBrghtG',
                   'GlobCluster', 'RadioG', 'EmissionG', 'BlueCompactG', 'StarburstG', 'HIIG', 'AGN',
                   'Seyfert', 'Seyfert1', 'Seyfert2', 'LINER', 'QSO', 'Blazar', 'BLLac']


def findClusters(saveToFile=False, name='clusterQuery.csv', numOfChildren=100):
    # Setting up query fields required

    customSimbad = Simbad()
    customSimbad.TIMEOUT = 1000
    customSimbad.add_votable_fields('z_value', 'ra', 'dec', 'dim_majaxis', 'dim_minaxis')
    customSimbad.remove_votable_fields('coordinates')
    qry = "children > " + str(numOfChildren) + " & otypes in ('ClG')"
    clusters = customSimbad.query_criteria(qry).to_pandas()
    clusters = clusters.dropna()


    if saveToFile:
        clusters.to_csv(name)

    return clusters


def getGalaxiesFromCluster(ra, dec, majaxis, clusterZ, name):
    global previousQueryTime

    clusterDir = "Data/GalaxySets/" + str(name) + ".csv"

    if not os.path.exists(clusterDir):
        # set up custom query and define which columns we want
        print("querying server")
        customSimbad = Simbad()
        customSimbad.TIMEOUT = 100

        customSimbad.add_votable_fields('ra(s)', 'dec(s)',
                                        'distance_result', 'otype', 'rv_value', 'z_value')
        customSimbad.remove_votable_fields('coordinates')

        my_radius = '0d' + str(int(majaxis)) + 'm0s'

        currentQueryTime = time.time()

        if currentQueryTime - previousQueryTime < 3:
            time.sleep(3 - currentQueryTime + previousQueryTime)

        try:
            table = customSimbad.query_region(coord.SkyCoord(ra, dec, unit=(u.hourangle, u.deg)),
                                              radius=my_radius)
        except:
            print("     ### Connection Timeout ###")
            time.sleep(30)
            table = None

        if table is None or len(table) == 0:
            return

        df = table.to_pandas()
        previousQueryTime = currentQueryTime
        df.to_csv(clusterDir)

    else:
        print("Found file")
        df = pd.read_csv(clusterDir)

    # Clean up data

    # This worked well for the a2029 galaxy
    # return gal_df[(gal_df['Z_VALUE'] < 0.09) * (gal_df['Z_VALUE'] > 0.065)]

    # Were going to cull the data, from our investigations into a2029, we found
    # z +- 0.013 is a good separation for the orbiting galaxies.

    redshiftDifferenceCutoff = 0.013

    gal_df = df[df['OTYPE'].isin(acceptableTypes)]
    gal_df = gal_df.dropna()
    gal_df = gal_df[np.abs(gal_df['Z_VALUE'] - clusterZ) <= redshiftDifferenceCutoff]

    # print(pd.Series({c: df[c].unique() for c in df})['OTYPE'])
    
    return gal_df

