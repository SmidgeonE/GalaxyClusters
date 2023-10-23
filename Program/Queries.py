from matplotlib import pyplot as plt
from Program.Equations import *
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
    currentQueryTime = time.time()

    if currentQueryTime - previousQueryTime < 3:
        print("waiting")
        time.sleep(3 - currentQueryTime + previousQueryTime)
        print("commencing")

    table = customSimbad.query_region(coord.SkyCoord(ra, dec, unit=(u.hourangle, u.deg)),
                                   radius=my_radius)

    if table is None or len(table) == 0:
        return

    df = table.to_pandas()
    previousQueryTime = currentQueryTime

    # Clean up data

    # This worked well for the a2029 galaxy
    # return gal_df[(gal_df['Z_VALUE'] < 0.09) * (gal_df['Z_VALUE'] > 0.065)]

    # Were going to cull the data, from our investigations into a2029, we found
    # z +- 0.013 is a good separation for the orbiting galaxies.

    redshiftDifferenceCutoff = 0.013

    print(pd.Series({c: df[c].unique() for c in df})['OTYPE'])

    # gal_df = df[df['OTYPE'] == 'Galaxy']
    gal_df = df

    gal_df = gal_df.dropna()
    gal_df = gal_df[np.abs(gal_df['Z_VALUE'] - clusterZ) <= redshiftDifferenceCutoff]

    return gal_df