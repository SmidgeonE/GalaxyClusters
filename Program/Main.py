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
    if makeGraphs:
        plt.show()

    # Now we plot a gaussian on the histogram for recession vel

    recessionVel = galaxiesOfCluster['RV_VALUE']
    fig, ax = plt.subplots(1, 1)
    histData = ax.hist(recessionVel, bins=numHistBins)

    velocityRange = np.linspace(np.min(recessionVel), np.max(recessionVel), len(recessionVel))
    gaussVals, sd = gauss(recessionVel, np.max(histData[0]), velocityRange)

    # Doing Chi squared values, we create an expected values list using the gauss function
    expectedValsRange = np.linspace(np.min(recessionVel), np.max(recessionVel), numHistBins)
    expectedVals = gauss(recessionVel, np.max(histData[0]), expectedValsRange)[0]

    ChiVal = np.sum(np.square(histData[0] - expectedVals) / expectedValsRange)
    print ("chi val : " + str(ChiVal))
    reducedChi = ChiVal / numHistBins

    ax.set_ylabel('num')
    ax.set_xlabel('RV_VALUE')
    ax.set_title('$R_v$ Values for ' + nameOfGalaxy + r"    $ \chi ^2$ = " + str(round(reducedChi, 2)))

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


## MAIN PROGRAM LOOP

# c = findClusters(saveToFile=True)

clustersSet = pd.read_csv('Data/clusterQuery.csv')

for i in clustersSet.index:
    galaxiesInCluster = getGalaxiesFromCluster(clustersSet['RA'][i],
                                               clustersSet['DEC'][i],
                                               clustersSet['GALDIM_MAJAXIS'][i],
                                               clustersSet['Z_VALUE'][i])

    if galaxiesInCluster is None or len(galaxiesInCluster.index) < 40:
        print("\nInsufficient Number of galaxies in cluster : " + str(clustersSet['MAIN_ID'][i]))
        continue

    CalculateAndPlotOmegaM(galaxiesInCluster, nameOfGalaxy=clustersSet['MAIN_ID'][i], makeGraphs=True)
    print("")


print("\nFinished Queries")
