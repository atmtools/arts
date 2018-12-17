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

N = 26

typ = ["s0", "f0", "SELF-G0-X0", "AIR-G0-X0", "SELF-G0-X1", "AIR-G0-X1", "SELF-D0-X0", "AIR-D0-X0", "SELF-D0-X1", "AIR-D0-X1", "SELF-G2-X0", "AIR-G2-X0", "SELF-G2-X1", "AIR-G2-X1", "SELF-D2-X0", "AIR-D2-X0", "SELF-D2-X1", "AIR-D2-X1", "SELF-FVC-X0", "AIR-FVC-X0", "SELF-FVC-X1", "AIR-FVC-X1", "SELF-ETA-X0", "AIR-ETA-X0", "SELF-ETA-X1", "AIR-ETA-X1"]

pert = [1e-17, 1e1, 20, 20, 0.8e-3, 0.8e-3, 1, 1, 0.8e-3, 0.8e-3, 2, 4, 1e-2, 1e-2, 1, 5e-1, 1e-2, 1e-2, 20, 20, 2, 2, 1e-4, 1e-4, 1e-4, 1e-4]

plot_shape = True
    
y = typhon.arts.xml.load("comparedata/y.xml")
dy = typhon.arts.xml.load("comparedata/dy.xml")
NF = 100
f = np.linspace(85, 115, NF)

if plot_shape:
    plt.figure(figsize=(5, 5))
    plt.plot(f, y)
    plt.show()

X = int(np.sqrt(N))
Y = N//X + 1
plt.figure(figsize=(8*Y, 8*X))
for i in range(N):
    yd = typhon.arts.xml.load("comparedata/y-d" + typ[i]+".xml")
    plt.subplot(X, Y, i + 1)
    plt.subplots_adjust(hspace=1, wspace=1)
    plt.plot(f, abs(yd-y)/pert[i])
    plt.semilogy(f, abs(dy[:, i]))
    plt.legend(('p', 'a'), fontsize=6)
    plt.title('d'+typ[i])
    plt.xlabel('Frequency [GHz]')
    plt.ylabel('dy/dx')
plt.show()


plt.figure(figsize=(8*Y, 8*X))
for i in range(N):
    yd = typhon.arts.xml.load("comparedata/y-d" + typ[i]+".xml")
    plt.subplot(X, Y, i + 1)
    plt.subplots_adjust(hspace=1, wspace=1)
    plt.plot(f, 100*(yd-y)/pert[i]/dy[:, i] - 100, 'k-*')
    plt.title('d'+typ[i])
    plt.xlabel('Frequency [GHz]')
    plt.ylabel('rel. dy/dx [%]')
plt.show()
