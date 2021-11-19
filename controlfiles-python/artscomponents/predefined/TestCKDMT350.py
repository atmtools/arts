# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import pyarts

SAVE = False
PLOT = False
CMPR = True

ws = pyarts.workspace.Workspace()

ws.stokes_dim = 1
ws.rtp_temperature = 250.
ws.rtp_pressure = 5e4
ws.rtp_vmr = np.array([1e-4])

ws.Touch(ws.jacobian_quantities)
f = np.linspace(1e-3, 20_001, 101)
ws.f_grid = f * 29979245800

vals = []
for spec in ["H2O-SelfContCKDMT350", "H2O-ForeignContCKDMT350"]:
    ws.abs_speciesSet(species=[spec])
    
    ws.propmat_clearsky_agenda_checked = 1
    ws.abs_xsec_agenda_checked = 1
    
    ws.propmat_clearskyInit()
    ws.propmat_clearskyAddPredefined()
    vals.append(1.0 *  ws.propmat_clearsky.value.data.data.flatten())
    
    if PLOT:
        plt.semilogy(f, ws.propmat_clearsky.value.data.data.flatten(), label=spec)

if PLOT:
    plt.legend()

res = pyarts.classes.ArrayOfVector(vals)

if SAVE:
    res.savexml("test_data_ckdmt350.xml")

if CMPR:
    cmpr = pyarts.classes.ArrayOfVector()
    cmpr.readxml("test_data_ckdmt350.xml")
    for i in range(len(cmpr)):
        assert np.isclose(cmpr[i].data, res[i].data).all()