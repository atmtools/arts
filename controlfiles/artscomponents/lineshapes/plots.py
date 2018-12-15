#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 22 12:14:03 2018

Plots testdata/dtest-"+ls+"/ files.

@author: larsson
"""

import typhon
import matplotlib.pyplot as plt
import numpy as np

lineshapes = ["doppler", "lorentz", "voigt", "htp-vp", "sdvp", "htp"]

num_dplots = {"doppler": 5,
              "lorentz": 31,
              "voigt": 31,
              "htp-vp": 25,
              "sdvp": 21,
              "htp": 29,}

type_dplots = {"doppler": ["T", "f", "vmr", "s0", "f0"],
               "lorentz": ["T", "f", "vmr", "s0", "f0", "SELF-G0-X0", "AIR-G0-X0", "SELF-G0-X1", "AIR-G0-X1", "SELF-D0-X0", "AIR-D0-X0", "SELF-D0-X1", "AIR-D0-X1", "SELF-Y-X0", "AIR-Y-X0", "SELF-Y-X1", "AIR-Y-X1", "SELF-Y-X2", "AIR-Y-X2", "SELF-G-X0", "AIR-G-X0", "SELF-G-X1", "AIR-G-X1", "SELF-G-X2", "AIR-G-X2", "SELF-DV-X0", "AIR-DV-X0", "SELF-DV-X1", "AIR-DV-X1", "SELF-DV-X2", "AIR-DV-X2"],
               "voigt": ["T", "f", "vmr", "s0", "f0", "SELF-G0-X0", "AIR-G0-X0", "SELF-G0-X1", "AIR-G0-X1", "SELF-D0-X0", "AIR-D0-X0", "SELF-D0-X1", "AIR-D0-X1", "SELF-Y-X0", "AIR-Y-X0", "SELF-Y-X1", "AIR-Y-X1", "SELF-Y-X2", "AIR-Y-X2", "SELF-G-X0", "AIR-G-X0", "SELF-G-X1", "AIR-G-X1", "SELF-G-X2", "AIR-G-X2", "SELF-DV-X0", "AIR-DV-X0", "SELF-DV-X1", "AIR-DV-X1", "SELF-DV-X2", "AIR-DV-X2"],
               "htp-vp": ["T", "f", "vmr", "s0", "f0", "SELF-G0-X0", "AIR-G0-X0", "SELF-G0-X1", "AIR-G0-X1", "SELF-D0-X0", "AIR-D0-X0", "SELF-D0-X1", "AIR-D0-X1", "SELF-Y-X0", "AIR-Y-X0", "SELF-Y-X1", "AIR-Y-X1", "SELF-Y-X2", "AIR-Y-X2", "SELF-G-X0", "AIR-G-X0", "SELF-G-X1", "AIR-G-X1", "SELF-G-X2", "AIR-G-X2"],
               "sdvp": ["T", "f", "vmr", "s0", "f0", "SELF-G0-X0", "AIR-G0-X0", "SELF-G0-X1", "AIR-G0-X1", "SELF-D0-X0", "AIR-D0-X0", "SELF-D0-X1", "AIR-D0-X1", "SELF-G2-X0", "AIR-G2-X0", "SELF-G2-X1", "AIR-G2-X1", "SELF-D2-X0", "AIR-D2-X0", "SELF-D2-X1", "AIR-D2-X1"],
               "htp": ["T", "f", "vmr", "s0", "f0", "SELF-G0-X0", "AIR-G0-X0", "SELF-G0-X1", "AIR-G0-X1", "SELF-D0-X0", "AIR-D0-X0", "SELF-D0-X1", "AIR-D0-X1", "SELF-G2-X0", "AIR-G2-X0", "SELF-G2-X1", "AIR-G2-X1", "SELF-D2-X0", "AIR-D2-X0", "SELF-D2-X1", "AIR-D2-X1", "SELF-FVC-X0", "AIR-FVC-X0", "SELF-FVC-X1", "AIR-FVC-X1", "SELF-ETA-X0", "AIR-ETA-X0", "SELF-ETA-X1", "AIR-ETA-X1"],}

pert_dplots = {"doppler": [0.0001, 100, 0.0001, 1e-30, 1e1],
               "lorentz": [0.0001, 100, 0.0001, 1e-30, 1e1, 20, 20, 0.8e-3, 0.8e-3, 1, 1, 0.8e-3, 0.8e-3, 1e-10, 1e-10, 1e-12, 1e-12, 0.8e-3, 0.8e-3, 1e-14, 1e-14, 1e-16, 1e-16, 0.8e-3, 0.8e-3, 1e-2, 1e-2, 1e-4, 1e-4, 0.8e-3, 0.8e-3],
               "voigt": [0.0001, 100, 0.0001, 1e-30, 1e1, 20, 20, 0.8e-3, 0.8e-3, 1, 1, 0.8e-3, 0.8e-3, 1e-10, 1e-10, 1e-12, 1e-12, 0.8e-3, 0.8e-3, 1e-14, 1e-14, 1e-16, 1e-16, 0.8e-3, 0.8e-3, 1e-2, 1e-2, 1e-4, 1e-4, 0.8e-3, 0.8e-3],
               "htp-vp": [0.0001, 100, 0.0001, 1e-30, 1e1, 20, 20, 0.8e-3, 0.8e-3, 1, 1, 0.8e-3, 0.8e-3, 1e-10, 1e-10, 1e-12, 1e-12, 0.8e-3, 0.8e-3, 1e-14, 1e-14, 1e-16, 1e-16, 0.8e-3, 0.8e-3],
               "sdvp": [0.0001, 100, 0.0001, 1e-30, 1e1, 20, 20, 0.8e-3, 0.8e-3, 1, 1, 0.8e-3, 0.8e-3, 2, 4, 1e-2, 1e-2, 1, 5e-1, 1e-2, 1e-2],
               "htp": [0.0001, 100, 0.0001, 1e-30, 1e1, 20, 20, 0.8e-3, 0.8e-3, 1, 1, 0.8e-3, 0.8e-3, 2, 4, 1e-2, 1e-2, 1, 5e-1, 1e-2, 1e-2, 20, 20, 2, 2, 1e-4, 1e-4, 1e-4, 1e-4],}

plot_shape = True

for ls in lineshapes:
    N = num_dplots[ls]
    typ = type_dplots[ls]
    pert = pert_dplots[ls]
    
    pm = typhon.arts.xml.load("testdata/test-" + ls + 
                              "/propmat.xml")[0].data[0, 0, :, 0]
    adpm = typhon.arts.xml.load("testdata/test-"+ls+"/dpropmat.xml")
    f = np.linspace(90, 110, 1001)
    
    if plot_shape:
        plt.figure(figsize=(5, 5))
        plt.plot(f, pm)
        plt.title(ls.upper())
        plt.show()
    
    X = int(np.sqrt(N))
    Y = N//X + 1
    plt.figure(figsize=(4*Y, 4*X))
    for i in range(N):
        pmd = typhon.arts.xml.load("testdata/test-" +
                                   ls + "/propmat-d" + 
                                   typ[i]+".xml")[0].data[0, 0, :, 0]
        plt.subplot(X, Y, i + 1)
        plt.plot(f, abs(pmd-pm)/pert[i])
        plt.semilogy(f, abs(adpm[i].data[0, 0, :, 0]))
        plt.legend(("analytical", "perturbed"),loc=8)
        plt.title('d'+typ[i])
    plt.tight_layout()
    plt.show()

#pmd =typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dT.xml")
#plt.subplot(6, 6, 1)
#plt.semilogy(f,abs(pmd[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 0.0001)
#plt.title("Temperature")

