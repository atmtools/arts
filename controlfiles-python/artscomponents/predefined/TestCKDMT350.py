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

##################### SILLY ARTS DEFINITIONS
ws.atmosphere_dim = 3
ws.p_grid = np.array([1000000000000000000000000000000000000000000000000000000])
ws.lat_grid = np.array([10000000000000000000000000000000000000000000000000000])
ws.lon_grid = np.array([10000000000000000000000000000000000000000000000000000])
ws.Touch(ws.select_abs_species)
ws.Touch(ws.predefined_model_data)
##################### /SILLY ARTS DEFINITIONS

f = np.linspace(1e-3, 20_001, 101)
ws.f_grid = f * 29979245800

vals = []
dvals = []
manual_dvals = []
specs = ["H2O-SelfContCKDMT350", "H2O-ForeignContCKDMT350"]
for spec in specs:
    
    ws.jacobianInit()
    ws.jacobianAddTemperature(g1=ws.p_grid, g2=ws.lat_grid, g3=ws.lon_grid)
    ws.jacobianAddWind(g1=ws.p_grid, g2=ws.lat_grid, g3=ws.lon_grid)
    ws.jacobianAddAbsSpecies(g1=ws.p_grid, g2=ws.lat_grid, g3=ws.lon_grid, species="H2O-161", for_species_tag=0)
    ws.jacobianAddAbsSpecies(g1=ws.p_grid, g2=ws.lat_grid, g3=ws.lon_grid, species="H2O-181", for_species_tag=0)
    ws.jacobianAddAbsSpecies(g1=ws.p_grid, g2=ws.lat_grid, g3=ws.lon_grid, species="H2O-SelfContCKDMT350", for_species_tag=1)
    ws.jacobianAddAbsSpecies(g1=ws.p_grid, g2=ws.lat_grid, g3=ws.lon_grid, species="H2O-ForeignContCKDMT350", for_species_tag=1)

    ws.abs_speciesSet(species=[spec])
    
    ws.propmat_clearsky_agenda_checked = 1
    
    ws.propmat_clearskyInit()
    ws.propmat_clearskyAddPredefined()
    vals.append(1.0 *  np.array(ws.propmat_clearsky.value.data).flatten())
    
    dvals.append([])
    for x in ws.dpropmat_clearsky_dx.value:
        dvals[-1].append(1.0 * np.array(x.data).flatten())
    
    manual_dvals.append([])
    ws.jacobianInit()
    ws.rtp_temperature.value += 0.1
    ws.propmat_clearskyInit()
    ws.propmat_clearskyAddPredefined()
    manual_dvals[-1].append(10.0 *  (np.array(ws.propmat_clearsky.value.data).flatten() - vals[-1]))
    ws.rtp_temperature.value -= 0.1
    
    ws.f_grid.value += 0.1
    ws.propmat_clearskyInit()
    ws.propmat_clearskyAddPredefined()
    manual_dvals[-1].append(10. *  (np.array(ws.propmat_clearsky.value.data).flatten() - vals[-1]))
    ws.f_grid.value -= 0.1
    
    ws.rtp_vmr.value += 1e-6
    ws.propmat_clearskyInit()
    ws.propmat_clearskyAddPredefined()
    manual_dvals[-1].append(1e6 *  (np.array(ws.propmat_clearsky.value.data).flatten() - vals[-1]))
    ws.rtp_vmr.value -= 1e-6
    
if PLOT:
    for i in range(len(specs)):
        spec = specs[i]
        val = vals[i]
        plt.semilogy(f, val, label=spec)
    plt.legend()
    plt.show()

res = pyarts.arts.ArrayOfVector(vals)
dres = pyarts.arts.ArrayOfArrayOfVector(dvals)
manual_dres = pyarts.arts.ArrayOfArrayOfVector(manual_dvals)

if SAVE:
    res.savexml("test_data_ckdmt350.xml")
    dres.savexml("test_data_ckdmt350.deriv.xml")
    

if CMPR:
    cmpr = pyarts.arts.ArrayOfVector()
    cmpr.readxml("test_data_ckdmt350.xml")
    for i in range(len(cmpr)):
        assert np.isclose(cmpr[i], res[i]).all()
        
    dcmpr = pyarts.arts.ArrayOfArrayOfVector()
    dcmpr.readxml("test_data_ckdmt350.deriv.xml")
    for i in range(len(dcmpr)):
        for j in range(len(dcmpr[i])):
            assert np.isclose(dcmpr[i][j], dres[i][j]).all()

# Note: manual_dres is not compared against but kept around.  The derivative of
# water VMR for the self-model has a factor x**2, where x is the VMR itself.
# The special derivative (for_species_tag=1) only cares for one of these xs so
# the manual derivative will be twice as large as it 'should' be.
