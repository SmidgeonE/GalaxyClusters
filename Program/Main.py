import numpy as np
from Program.Queries import *


def CalculateAndPlotOmegaM(galaxiesOfCluster, nameOfGalaxy, makeGraphs=False):
    global histData
    global expectedVals
    numHistBins = 40

    # Now we generate some graphs of the cluster

    fig2, ax2 = plt.subplots(1, 1)
    ax2.hist(galaxiesOfCluster['Z_VALUE'], bins=numHistBins)
    ax2.set_ylabel('num')
    ax2.set_xlabel('Z_VALUE')
    ax2.set_title('Z Values for ' + nameOfGalaxy)

    # if makeGraphs:
    #     plt.show()

    # Now we plot a gaussian on the histogram for recession vel

    recessionVel = galaxiesOfCluster['RV_VALUE']
    fig, ax = plt.subplots(1, 2)
    histData = ax[0].hist(recessionVel, bins=numHistBins)
    ax[0].set_ylabel('num')
    ax[0].set_xlabel('RV_VALUE')

    velocityRange = np.linspace(np.min(recessionVel), np.max(recessionVel), len(recessionVel))
    gaussVals, sd = gauss(recessionVel, np.max(histData[0]), velocityRange)

    ax[0].plot(velocityRange, gaussVals)

    # Doing Chi squared values, we create an expected values list using the gauss function

    expectedValsRange = np.linspace(np.min(recessionVel), np.max(recessionVel), numHistBins)
    expectedVals = gauss(recessionVel, np.max(histData[0]), expectedValsRange)[0]

    # We need to cull data bins that have a low expectation value, as the chi value will be massive for even a small
    # deviation

    expectedVals = np.where(expectedVals < 1, np.nan, expectedVals)
    excludedDataMask = np.where(np.isnan(expectedVals), 1, 0)
    correctDataMask = np.where(np.isnan(expectedVals), 0, 1)
    excludedData = histData[0] * excludedDataMask
    correctedData = histData[0] * correctDataMask

    # Plotting the new culled data on top of the unculled data...

    edges = histData[1]
    centres = (edges[:-1] + edges[1:]) / 2
    barWidth = edges[1]-edges[0]
    ax[1].bar(centres, excludedData, width=barWidth, color='green')
    ax[1].set_title('Culled Data')
    ax[1].set_xlabel('RV_VALUE')
    ax[1].set_ylim([0, ax[0].get_ylim()[1]])

    # Calculating chi values, plotting it in the title

    ChiVal = np.nansum(np.square(correctedData - expectedVals) / expectedVals)
    reducedChi = ChiVal / numHistBins

    ax[0].set_title('$R_v$ Values for ' + nameOfGalaxy + r"    $ \chi ^2$ = " + str(round(reducedChi, 2)))

    if makeGraphs:
        plt.show()

    # Removing culled data from the data pool for calculations...
    numOfCulledGals = 0

    for point in excludedData:
        if point == 0:
            continue

        for index, galaxy in galaxiesOfCluster.iterrows():
            if np.abs(float(galaxy['RV_VALUE'])-point) <= barWidth / 2:
                galaxiesOfCluster = galaxiesOfCluster[galaxiesOfCluster.index != index]
                numOfCulledGals += 1


    # Recalculate standard deviation with new dataset...

    sd = np.nanstd(galaxiesOfCluster['RV_VALUE'])

    # Apply culled data to equations to find mass of cluster + omega_m

    averageZ = galaxiesOfCluster['Z_VALUE'].mean()
    M_200val = M_200(sd * 1000, averageZ)
    f_ourGas = M_gas(M_200val) / M_200val
    f_b = f_star + f_ourGas
    Omega_m = Omega_b / f_b

    print("\nDATA FOR : " + nameOfGalaxy)
    print("--- omega_rv: " + str(sd * 1000))
    print("--- M_200 value: " + str(M_200val))
    print("--- omega m : " + str(Omega_m))

    return Omega_m


## MAIN PROGRAM LOOP

# c = findClusters(saveToFile=True)

clustersSet = pd.read_csv('Data/clusterQuery.csv')

for i in clustersSet.index:
    # if clustersSet['MAIN_ID'][i] != "M  60":
    #     continue


    galaxiesInCluster = getGalaxiesFromCluster(clustersSet['RA'][i],
                                               clustersSet['DEC'][i],
                                               clustersSet['GALDIM_MAJAXIS'][i],
                                               clustersSet['Z_VALUE'][i])

    if galaxiesInCluster is None:
        print("\nCluster Data is empty!")
        continue

    if len(galaxiesInCluster.index) < 20:
        print("\nInsufficient Number of galaxies in cluster : " + str(clustersSet['MAIN_ID'][i] +
                                                                      ". \nNum of galaxies in cluster:"))
        print(str(len(galaxiesInCluster.index)))
        continue

    CalculateAndPlotOmegaM(galaxiesInCluster, nameOfGalaxy=clustersSet['MAIN_ID'][i], makeGraphs=True)
    print("------------------------\n")

    break


print("\nFinished Queries")
