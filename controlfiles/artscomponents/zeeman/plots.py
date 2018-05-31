#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 22 12:14:03 2018

Plots testdata/dtest-"+ls+"/ files.

First figure shows values of analytically and perturbed
computations of the derivatives for a few variable.

Second figure shows the relative differences between
analytically and perturbed calculations for the same
derivatives.

@author: larsson
"""

import typhon
import numpy as np
import matplotlib.pyplot as plt

###
# SELECT PLOT TYPE BY NAMING ls as one of below
###


f = np.linspace(-1, 1, 501)
pm = typhon.arts.xml.load("testdata/zeeman/propmat.xml")[0].data[0, 0]

plt.figure(figsize=(7, 7))
plt.subplot(4, 4, 1)
plt.plot(f, pm[:, 0])
plt.subplot(4, 4, 2)
plt.plot(f, pm[:, 1])
plt.subplot(4, 4, 3)
plt.plot(f, pm[:, 2])
plt.subplot(4, 4, 4)
plt.plot(f, pm[:, 3])
plt.subplot(4, 4, 5)
plt.plot(f, pm[:, 1])
plt.subplot(4, 4, 6)
plt.plot(f, pm[:, 0])
plt.subplot(4, 4, 7)
plt.plot(f, pm[:, 4])
plt.subplot(4, 4, 8)
plt.plot(f, pm[:, 5])
plt.subplot(4, 4, 9)
plt.plot(f, pm[:, 2])
plt.subplot(4, 4, 10)
plt.plot(f, -pm[:, 4])
plt.subplot(4, 4, 11)
plt.plot(f, pm[:, 0])
plt.subplot(4, 4, 12)
plt.plot(f, pm[:, 6])
plt.subplot(4, 4, 13)
plt.plot(f, pm[:, 3])
plt.subplot(4, 4, 14)
plt.plot(f, -pm[:, 5])
plt.subplot(4, 4, 15)
plt.plot(f, -pm[:, 6])
plt.subplot(4, 4, 16)
plt.plot(f, pm[:, 0])
plt.tight_layout()
