#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 12:35:24 2023

@author: joshsaunders
"""

import numpy as np

repeats = 1000
sigmaV = np.zeros(repeats)


def BootstrapErr(cluster):
    for repeat in range(repeats):
        recessionVel = cluster['RV_VALUE']

        # get the number of velocities in the array
        N = len(recessionVel)

        # sample the same number of velocities from the original array with replacement
        recessionVel_sample = np.random.choice(recessionVel, N, replace=True)

        sigmaV[repeat] = np.std(recessionVel_sample)

    return np.std(sigmaV)
