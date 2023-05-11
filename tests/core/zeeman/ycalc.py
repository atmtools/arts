#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 09:49:38 2021

@author: larsson
"""

import pyarts
import os
import numpy as np

arts = pyarts.workspace.Workspace()

# %% Options

ATMPATH = "planets/Earth/Fascod/tropical2/tropical"
LINEPATH = "line"
DYNMAG = True

NF = 100
NR = 9

SHOW_PLOTS = True

DO_SAVE = False

CF_SAVE = True

# %% Optional imports

if SHOW_PLOTS:
    import matplotlib.pyplot as plt

# %% Optional values

MAGSTR = 45e-6  # T
MAGTHE = 30  # deg
MAGETA = 45  # deg
CENTRAL_LINE_FREQ = 118750348044.712
SPECTROMETER_HW = 5e6

# %% Silly ARTS-isms
    
arts.abs_f_interp_order = 100_000_000
arts.Touch(arts.iy_aux_vars)
arts.Touch(arts.surface_props_names)
arts.Touch(arts.transmitter_pos)

# %% Agendas
arts.iy_main_agendaSetByPart(rte_option="Emission",
                             ppath_option="Geometric",
                             background_option="Emission",
                             propagation_properties_option="FromPropmat",
                             )

@pyarts.workspace.arts_agenda(ws=arts, set_agenda=True)
def surface_rtprop_agenda(arts):
    arts.InterpSurfaceFieldToPosition()
    arts.surfaceBlackbody()

@pyarts.workspace.arts_agenda(ws=arts)
def iy_space_agenda_cosmic_background(arts):
    arts.Ignore(arts.rtp_pos)
    arts.Ignore(arts.rtp_los)
    arts.MatrixCBR(arts.iy, arts.stokes_dim, arts.f_grid)
arts.Copy(arts.iy_space_agenda, iy_space_agenda_cosmic_background)

@pyarts.workspace.arts_agenda(ws=arts)
def iy_surface_agenda(arts):
    arts.SurfaceDummy()
    arts.InterpSurfaceFieldToPosition()
    arts.iySurfaceRtpropAgenda()
arts.Copy(arts.iy_surface_agenda, iy_surface_agenda)

@pyarts.workspace.arts_agenda(ws=arts)
def water_psat_agenda(arts):
    arts.water_p_eq_fieldMK05()
arts.Copy(arts.water_p_eq_agenda, water_psat_agenda)
    
# %% Calculations

arts.jacobianOff()
arts.iy_unit = "PlanckBT"
arts.ppath_lstep = 1000.
arts.ppath_lmax = 10e3
arts.ppath_lraytrace = 1e2
arts.rte_alonglos_v = 0.0
arts.nlteOff()

# %% Species and line absorption

arts.abs_speciesSet(species=[f"O2-Z-66-{CENTRAL_LINE_FREQ-1}-{CENTRAL_LINE_FREQ+1}"])
arts.abs_lines_per_speciesReadSpeciesSplitCatalog(basename = "lines/")
arts.Wigner6Init()

# %% Use the automatic agenda setter
arts.propmat_clearsky_agendaAuto(manual_mag_field=not DYNMAG,
                                         H=MAGSTR,
                                         theta=MAGTHE,
                                         eta=MAGETA)

# %% Grids and planet

arts.refellipsoidEarth(model = "Sphere")
arts.surface_fieldInit()
arts.surface_fieldSet(value=0.0, key="h")
arts.surface_fieldSet(value=200.0, key="t")

# %% Atmosphere
arts.Touch(arts.time)
arts.atm_fieldInit(toa=95e3)
arts.atm_fieldAddNumericData(key=arts.abs_species.value[0], data=0.21)
arts.atm_fieldAddCustomDataFile(key="t", filename=ATMPATH+".t.xml")
arts.atm_fieldAddCustomDataFile(key="p", filename=ATMPATH+".p.xml")
arts.atm_fieldIGRF()

# %% Sensor

arts.stokes_dim = 4
arts.f_grid = np.linspace(-SPECTROMETER_HW, SPECTROMETER_HW, NF) + CENTRAL_LINE_FREQ
arts.sensor_pos = np.zeros((NR, 3))
arts.sensor_pos.value.value[:, 0] = 300e3
arts.sensor_pos.value.value[:, 1] = np.linspace(-80, 80, NR)
arts.sensor_los = np.zeros((NR, 2))
arts.sensor_los.value.value[:, 0] = 180
arts.sensorOff()

# %% Computechecks

arts.atmgeom_checked = 1
arts.lbl_checkedCalc()
arts.atmfields_checked = 1
arts.cloudbox_checked = 1
arts.cloudbox_on = 0
arts.sensor_checkedCalc()
arts.propmat_clearsky_agenda_checkedCalc()

# %% SURFACE
arts.surface_elevation = pyarts.arts.GriddedField2([[0], [0]], [[0]], ["Latitude", "Longitude"])
arts.refellipsoidEarthZZZ(refellipsoidZZZ=arts.refellipsoid, model="Sphere")

arts.Touch(arts.iy_cloudbox_agenda)

# %% Compute 
print(arts.iy_main_agenda.value)

pyarts.arts.omp_set_num_threads(1)
arts.yCalc()

arts.rte_los = arts.sensor_los.value[0]
arts.rte_pos = arts.sensor_pos.value[0]
arts.ppathGeometric()
arts.ppvar_atmFromPath()
atm = arts.ppvar_atm.value
path = arts.ppath.value

plt.plot(arts.y.value.flatten()[::4])
arts.VectorCreate("yref")
arts.ReadXML(arts.yref, "refdata.xml")
plt.plot(arts.yref.value.flatten()[::4])
plt.show()
