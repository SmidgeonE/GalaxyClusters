import numpy as np

from Program.Constants import Omega_b, f_star
from Program.Queries import *
from Program.bootstrap import *
import matplotlib.pyplot as plt
import pandas as pd
import Program.MonteCarlo as mc
from Program.bootstrap import BootstrapErr
import scienceplots

plt.style.use(['science', 'notebook', 'grid'])


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
    sigmaRvError = BootstrapErr(galaxiesOfCluster)

    M_200val = M_200(sd * 1000, averageZ)
    M_200err = mc.M_200Err(averageZ, sd * 1000, sigmaRvError * 1000)

    f_ourGas = F_gas(averageZ, M_200val)
    f_ourGasErr = mc.f_gasErr(averageZ, M_200val, M_200err)

    f_b = f_star + f_ourGas
    Omega_m = Omega_b / f_b


    print("\nDATA FOR: " + nameOfGalaxy)
    print("--- f_gas: " + str(f_ourGas))
    print("--- f_gas error:" + str(f_ourGasErr))
    print("--- rv mean: " + str(np.mean(recessionVel)))
    print("--- sigma_rv: " + str(sd))
    print("--- sigma_rv error: " + str(sigmaRvError))
    print("--- M_200 value: " + str(M_200val))
    print("--- M_200 error: " + str(M_200err))
    print("--- Omega_M: " + str(Omega_m))
    print("--- Omega_M MC error: " + str(mc.OmegaMErrorGeneral(f_ourGas, f_ourGasErr)))

    # If the Chi value is too large, we will ignore this data set

    if reducedChi > 0.8:
        return np.nan, np.nan, np.nan, np.nan

    return Omega_m, mc.OmegaMErrorGeneral(f_ourGas, f_ourGasErr), M_200val, M_200err


def RunAnalysisForSet(clusterDir, returnRawVals=False):
    clustersSet = pd.read_csv(clusterDir)
    omega_M = np.zeros(len(clustersSet.index))
    omega_Merrs = np.zeros(len(clustersSet.index))
    M_200 = np.zeros(len(clustersSet.index))
    M_200err = np.zeros(len(clustersSet.index))

    for i in clustersSet.index:
        # if clustersSet['MAIN_ID'][i] != "ACO  2029":
        #     continue
        print("currently reading cluster " + str(i))

        galaxiesInCluster = getGalaxiesFromCluster(clustersSet['RA'][i],
                                                   clustersSet['DEC'][i],
                                                   clustersSet['GALDIM_MAJAXIS'][i],
                                                   clustersSet['Z_VALUE'][i],
                                                   clustersSet['MAIN_ID'][i])

        if galaxiesInCluster is None:
            print("\nCluster Data is empty!")
            print("------------------------\n")
            continue

        if len(galaxiesInCluster.index) < 20:
            print("\nInsufficient Number of galaxies in cluster : " + str(clustersSet['MAIN_ID'][i] +
                                                                          ". \nNum of galaxies in cluster:"))
            print(str(len(galaxiesInCluster.index)))
            print("------------------------\n")
            continue


        omega_M[i], omega_Merrs[i], M_200[i], M_200err[i] = CalculateAndPlotOmegaM(galaxiesInCluster,
                                            nameOfGalaxy=clustersSet['MAIN_ID'][i],
                                            makeGraphs=False)
        print("------------------------\n")

    # If we skipped over some values, we need to ensure we're not counting the leftover zeros
    omega_M[omega_M == 0] = np.nan

    print("\nFinished Queries")
    print("Final Results: ")
    print("--- Overall Omega_M: " + str(np.nanmean(omega_M)))
    print("--- Overall Omega_M standard deviation: " + str(np.nanstd(omega_M)))
    print("--- Number of Clusters Analysed: " + str(len(clustersSet.index)))
    print("--- Number of Clusters Accepted: " + str(sum(~np.isnan(omega_M))))

    if returnRawVals:
        return omega_M, omega_Merrs, M_200, M_200err

    return np.nanmean(omega_M), np.nanstd(omega_M), \
           len(clustersSet.index), sum(~np.isnan(omega_M))




## MAIN PROGRAM LOOP


dirs = ["Data/clusterQuery0children.csv",
        "Data/clusterQuery1children.csv",
        "Data/clusterQuery20children.csv",
        "Data/clusterQuery21children.csv",
        "Data/clusterQuery40children.csv",
        "Data/clusterQuery50children.csv",
        "Data/clusterQuery100children.csv"]


def PlotOmegaMAgainstChiAndChildren():
    # This was using cutoff chi = 0.8

    # resultsForDataSets = np.zeros((len(dirs), 4))
    #
    # for i, d in enumerate(dirs):
    #     resultsForDataSets[i, :] = RunAnalysisForSet(d)

    dfBetterChi = pd.read_csv("Data/Results/OmegaMForChi0.8.csv")
    dfWorseChi = pd.read_csv("Data/Results/OmegaMForChi1.6.csv")

    # Changing to string array for better axis spacing
    strings = []
    for item in dfBetterChi["numChildren"].values:
        strings.append(str(item))

    plt.errorbar(strings, dfBetterChi["OmegaM"].values, marker='x', capsize=5, yerr=dfBetterChi['sd'], alpha=0.7)
    plt.errorbar(strings, dfWorseChi["OmegaM"].values, marker='x', capsize=5, yerr=dfWorseChi['sd'], alpha=0.7)
    plt.ylabel(r'$\Omega_{M}$')
    plt.xlabel('Number of Children Required')
    plt.legend([r"$\chi < 0.8$", r"$\chi < 1.6$"])
    plt.show()


def PlotOmegaMAgainstMass():
    data = pd.read_csv("Data/Results/OmegaMAgainstMtot.csv").transpose()
    plt.errorbar(data[0][1:], data[2][1:] / 10E46, marker='x', capsize=5, yerr=data[3][1:] / 10E46, xerr=data[1][1:], alpha=0.7)
    plt.ylabel(r'$M_{tot} / kg / 10^{46}$')
    plt.xlabel(r'$\Omega_{M}$')
    plt.legend([r"$\chi < 0.8$"])
    plt.show()

PlotOmegaMAgainstMass()


