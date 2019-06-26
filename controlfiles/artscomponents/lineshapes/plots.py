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

SAVE_FIGURES = False

lineshapes = ["lp", "lm-lp", "vp", "lm-vp",
              "htp-vp", 'htp-sdvp', 'htp']

freq = np.linspace(90, 110, 101) - 100

derivs = ['dT', 'df', 'dvmr', 'ds0', 'df0',
          'G0-X0', 'G0-X1', 'D0-X0', 'D0-X1', 'G2-X0', 'G2-X1', 'D2-X0', 'D2-X1',
          'FVC-X0', 'FVC-X1', 'ETA-X0', 'ETA-X1', 'Y-X0', 'Y-X1', 'Y-X2',
          'G-X0', 'G-X1', 'G-X2', 'DV-X0', 'DV-X1', 'DV-X2']
perturbs = [0.0001, 100, 0.0001, 1e-17, 1e1, 20, 0.8e-3, 1, 0.8e-3, 4, 1e-2, 5e-1, 1e-2, 20, 2, 1e-6, 1e-8,
            1e-10, 1e-12, 0.8e-3, 1e-14, 1e-16, 0.8e-3, 1e-2, 1e-4, 0.8e-3]

nshapes = len(lineshapes)
nderivs = len(derivs)

# Plot the lineshape
plt.figure(figsize=(8, 8))
plt.clf()
for il in range(nshapes):
    ls = lineshapes[il]
    x = typhon.arts.xml.load('testdata/test-{}/propmat.xml'.format(ls))[0].data.flatten()
    plt.semilogy(freq, x)
    plt.legend(lineshapes)
    plt.xlabel('Freq [GHZ]')
    plt.ylabel('Absorption Single Line [m$^{-1}$]')
plt.title('High pressure lineshapes')
if SAVE_FIGURES:
    plt.savefig('ls.png', dpi=300)
else:
    plt.show()


# Plot the lineshape diff to VP
plt.figure(figsize=(8, 8))
plt.clf()
for il in range(nshapes):
    ls = lineshapes[il]
    x = typhon.arts.xml.load('testdata/test-{}/propmat.xml'.format(ls))[0].data.flatten()
    x0 = typhon.arts.xml.load('testdata/test-{}/propmat.xml'.format('vp'))[0].data.flatten()
    plt.plot(freq, x-x0)
    plt.legend(lineshapes)
    plt.xlabel('Freq [GHZ]')
    plt.ylabel('Absorption Single Line Diff [m$^{-1}$]')
plt.title('Diff of high pressure lineshapes to Voigt')
if SAVE_FIGURES:
    plt.savefig('ls-diff.png', dpi=300)
else:
    plt.show()


# Plot the lineshape
X = int(np.sqrt(nderivs))
Y = nderivs//X + 1
plt.figure(figsize=(10*X, 8*Y))
plt.clf()
dpmat = np.zeros((nshapes), dtype=int)
for ids in range(nderivs):
    plt.subplot(X, Y, ids+1)
    dx = perturbs[ids]
    der = derivs[ids]
    this_leg = []
    for il in range(nshapes):
        ls = lineshapes[il]
        x0 = typhon.arts.xml.load('testdata/test-{}/propmat.xml'.format(ls))[0].data.flatten()
        try:
            x1 = typhon.arts.xml.load('testdata/test-{}/propmat-{}.xml'.format(ls, der))[0].data.flatten()
        except FileNotFoundError:
            continue
        dp = typhon.arts.xml.load('testdata/test-{}/dpropmat.xml'.format(ls))[dpmat[il]].data.flatten()
        plt.plot(freq, (x1-x0) / dx)
        plt.plot(freq, dp)
        plt.xlabel('Freq [GHZ]')
        this_leg.append('{} perturbed'.format(ls))
        this_leg.append('{} analytical'.format(ls))
        dpmat[il] += 1
    plt.legend(this_leg)
    plt.title(der)
plt.tight_layout()
if SAVE_FIGURES:
    plt.savefig('dls.png', dpi=100)
else:
    plt.show()


# Plot the lineshape
X = int(np.sqrt(nderivs))
Y = nderivs//X + 1
plt.figure(figsize=(10*X, 8*Y))
plt.clf()
dpmat = np.zeros((nshapes), dtype=int)
for ids in range(nderivs):
    plt.subplot(X, Y, ids+1)
    dx = perturbs[ids]
    der = derivs[ids]
    this_leg = []
    for il in range(nshapes):
        ls = lineshapes[il]
        x0 = typhon.arts.xml.load('testdata/test-{}/propmat.xml'.format(ls))[0].data.flatten()
        try:
            x1 = typhon.arts.xml.load('testdata/test-{}/propmat-{}.xml'.format(ls, der))[0].data.flatten()
        except FileNotFoundError:
            continue
        dp = typhon.arts.xml.load('testdata/test-{}/dpropmat.xml'.format(ls))[dpmat[il]].data.flatten()
        plt.plot(freq, (x1-x0) / dx / dp)
        plt.xlabel('Freq [GHZ]')
        this_leg.append('{} ratio'.format(ls))
        dpmat[il] += 1
    plt.legend(this_leg)
    plt.title(der)
plt.tight_layout()
if SAVE_FIGURES:
    plt.savefig('dls-ratio.png', dpi=100)
else:
    plt.show()