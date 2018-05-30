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

ls = "doppler"
ls = "lorentz"
ls = "fake-htp"
ls = "voigt"


pm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat.xml")
adpm = typhon.arts.xml.load("testdata/test-"+ls+"/dpropmat.xml")
f = np.linspace(90, 110, 1001)

labels = ["Propmat", "T", "Freq", "VMR", "S0", "F0",
          "PB-SG", "PB-FG", "PB-SE", "PB-FE", "PB-DV",
          "LM-Y0", "LM-G0", "LM-DV0",
          "LM-Y1", "LM-G1", "LM-DV1",
          "LM-YE", "LM-GE", "LM-DVE"]

plt.figure(figsize=(15, 8))
plt.subplot(4, 5, 1)
plt.title(labels[0])
plt.plot(f, pm[0].data[0, 0, :, 0])
i = 2
for ipm in adpm:
    plt.subplot(4, 5, i)
    plt.title("d" + labels[i-1])
    plt.plot(f, ipm.data[0, 0, :, 0])
    i += 1
plt.tight_layout()

plt.subplot(4, 5, 1)
plt.plot(f, pm[0].data[0, 0, :, 0])
i = 2

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dT.xml")
plt.subplot(4, 5, i)
plt.title("d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 0.0001)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-df.xml")
plt.subplot(4, 5, i)
plt.title("d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e2)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dvmr.xml")
plt.subplot(4, 5, i)
plt.title("d" + labels[i-1])
i += 1
# Note: difference here is if VMR = 20%, then ~20% of derivative is lineshape,
# so we cannot retrieve number densities for non-trace gasses before we
# implement lineshape derivatives for VMR...
# This holds true for all retrievals presently in ARTS...
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 0.0001)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-ds0.xml")
plt.subplot(4, 5, i)
plt.title("d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-30)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-df0.xml")
plt.subplot(4, 5, i)
plt.title("d" + labels[i-1])
i += 1
# Note: the cause of the slope is unknown.  Its amplitude depends on how
# accurate df0 is though, so the effect is minor
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e1)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dpb-sg.xml")
plt.subplot(4, 5, i)
plt.title("d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e1)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dpb-fg.xml")
plt.subplot(4, 5, i)
plt.title("d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e1)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dpb-se.xml")
plt.subplot(4, 5, i)
plt.title("d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-5)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dpb-fe.xml")
plt.subplot(4, 5, i)
plt.title("d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-5)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dpb-fs.xml")
plt.subplot(4, 5, i)
plt.title("d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e0)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlm-y0.xml")
plt.subplot(4, 5, i)
plt.title("d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-10)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlm-g0.xml")
plt.subplot(4, 5, i)
plt.title("d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-14)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlm-df0.xml")
plt.subplot(4, 5, i)
plt.title("d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-3)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlm-y1.xml")
plt.subplot(4, 5, i)
plt.title("d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-12)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlm-g1.xml")
plt.subplot(4, 5, i)
plt.title("d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-16)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlm-df1.xml")
plt.subplot(4, 5, i)
plt.title("d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-5)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlm-ye.xml")
plt.subplot(4, 5, i)
plt.title("d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-3)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlm-ge.xml")
plt.subplot(4, 5, i)
plt.title("d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-3)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlm-dfe.xml")
plt.subplot(4, 5, i)
plt.title("d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-3)

plt.tight_layout()
plt.show()

plt.figure(figsize=(15, 8))
dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dT.xml")
i=2
plt.subplot(4, 5, i)
plt.title("rel. d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 0.0001 / adpm[i-3].data[0,0,:,0] - 1)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-df.xml")
plt.subplot(4, 5, i)
plt.title("rel. d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e2 / adpm[i-3].data[0,0,:,0] - 1)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dvmr.xml")
plt.subplot(4, 5, i)
plt.title("rel. d" + labels[i-1])
i += 1
# Note: difference here is if VMR = 20%, then ~20% of derivative is lineshape,
# so we cannot retrieve number densities for non-trace gasses before we
# implement lineshape derivatives for VMR...
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 0.0001 / adpm[i-3].data[0,0,:,0] - 1)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-ds0.xml")
plt.subplot(4, 5, i)
plt.title("rel. d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-30 / adpm[i-3].data[0,0,:,0] - 1)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-df0.xml")
plt.subplot(4, 5, i)
plt.title("rel. d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e1 / adpm[i-3].data[0,0,:,0] - 1)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dpb-sg.xml")
plt.subplot(4, 5, i)
plt.title("rel. d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e1 / adpm[i-3].data[0,0,:,0] - 1)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dpb-fg.xml")
plt.subplot(4, 5, i)
plt.title("rel. d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e1 / adpm[i-3].data[0,0,:,0] - 1)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dpb-se.xml")
plt.subplot(4, 5, i)
plt.title("rel. d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-5 / adpm[i-3].data[0,0,:,0] - 1)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dpb-fe.xml")
plt.subplot(4, 5, i)
plt.title("rel. d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-5 / adpm[i-3].data[0,0,:,0] - 1)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dpb-fs.xml")
plt.subplot(4, 5, i)
plt.title("rel. d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e0 / adpm[i-3].data[0,0,:,0] - 1)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlm-y0.xml")
plt.subplot(4, 5, i)
plt.title("rel. d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-10 / adpm[i-3].data[0,0,:,0] - 1)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlm-g0.xml")
plt.subplot(4, 5, i)
plt.title("rel. d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-14 / adpm[i-3].data[0,0,:,0] - 1)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlm-df0.xml")
plt.subplot(4, 5, i)
plt.title("rel. d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-3 / adpm[i-3].data[0,0,:,0] - 1)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlm-y1.xml")
plt.subplot(4, 5, i)
plt.title("rel. d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-12 / adpm[i-3].data[0,0,:,0] - 1)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlm-g1.xml")
plt.subplot(4, 5, i)
plt.title("rel. d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-16 / adpm[i-3].data[0,0,:,0] - 1)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlm-df1.xml")
plt.subplot(4, 5, i)
plt.title("rel. d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-5 / adpm[i-3].data[0,0,:,0] - 1)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlm-ye.xml")
plt.subplot(4, 5, i)
plt.title("rel. d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-3 / adpm[i-3].data[0,0,:,0] - 1)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlm-ge.xml")
plt.subplot(4, 5, i)
plt.title("rel. d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-3 / adpm[i-3].data[0,0,:,0] - 1)

dpm = typhon.arts.xml.load("testdata/test-"+ls+"/propmat-dlm-dfe.xml")
plt.subplot(4, 5, i)
plt.title("rel. d" + labels[i-1])
i += 1
plt.plot(f, (dpm[0].data[0, 0, :, 0] - pm[0].data[0, 0, :, 0]) / 1e-3 / adpm[i-3].data[0,0,:,0] - 1)

plt.tight_layout()
plt.show()