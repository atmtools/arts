#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 09:49:38 2021

@author: larsson
"""

import pyarts
import os
import numpy as np

testdir = os.path.dirname(__file__)

arts = pyarts.workspace.Workspace()
arts.verbosityInit()

# %% Options

ATMPATH = "atm"
LINEPATH = "line"
DYNMAG = True

NF = 100
NR = 9

SHOW_PLOTS = False

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
    
arts.abs_f_interp_order = 100_000_000_000_000_000_000_000_000_000_000_000_000
arts.Touch(arts.iy_aux_vars)
arts.Touch(arts.surface_props_data)
arts.Touch(arts.surface_props_names)
arts.Touch(arts.transmitter_pos)

# %% Agendas

@pyarts.workspace.arts_agenda
def propmat_clearsky_agenda_dyn(ws):
    ws.propmat_clearskyInit()
    ws.propmat_clearskyAddZeeman()

@pyarts.workspace.arts_agenda
def propmat_clearsky_agenda_nodyn(ws):
    ws.propmat_clearskyInit()
    ws.propmat_clearskyAddZeeman(manual_zeeman_tag=True,
                                 manual_zeeman_magnetic_field_strength=MAGSTR,
                                 manual_zeeman_theta=MAGTHE,
                                 manual_zeeman_eta=MAGETA)
arts.Copy(arts.propmat_clearsky_agenda,
          propmat_clearsky_agenda_dyn if DYNMAG else propmat_clearsky_agenda_nodyn)

@pyarts.workspace.arts_agenda
def ppath_agenda_step_by_step(ws):
    ws.Ignore(ws.rte_pos2)
    ws.ppathStepByStep()
arts.Copy(arts.ppath_agenda, ppath_agenda_step_by_step)

@pyarts.workspace.arts_agenda
def iy_main_agenda_emission(ws):
    ws.ppathCalc()
    ws.iyEmissionStandard()
arts.Copy(arts.iy_main_agenda, iy_main_agenda_emission)

@pyarts.workspace.arts_agenda
def surface_rtprop_agenda(ws):
    ws.InterpSurfaceFieldToPosition(out=ws.surface_skin_t, field=ws.t_surface)
    ws.surfaceBlackbody()
arts.Copy(arts.surface_rtprop_agenda, surface_rtprop_agenda)
 
@pyarts.workspace.arts_agenda
def ppath_step_agenda_geometric(ws):
    ws.Ignore(ws.t_field)
    ws.Ignore(ws.vmr_field)
    ws.Ignore(ws.f_grid)
    ws.Ignore(ws.ppath_lraytrace)
    ws.ppath_stepGeometric()
arts.Copy(arts.ppath_step_agenda, ppath_step_agenda_geometric)

@pyarts.workspace.arts_agenda
def iy_space_agenda_cosmic_background(ws):
    ws.Ignore(ws.rtp_pos)
    ws.Ignore(ws.rtp_los)
    ws.MatrixCBR(ws.iy, ws.stokes_dim, ws.f_grid)
arts.Copy(arts.iy_space_agenda, iy_space_agenda_cosmic_background)

@pyarts.workspace.arts_agenda
def geo_pos_agenda(ws):
    ws.Ignore(ws.ppath)
    ws.VectorSet(ws.geo_pos, np.array([]))
arts.Copy(arts.geo_pos_agenda, geo_pos_agenda)

@pyarts.workspace.arts_agenda
def iy_surface_agenda(ws):
    ws.SurfaceDummy()
    ws.iySurfaceRtpropAgenda()
arts.Copy(arts.iy_surface_agenda, iy_surface_agenda)

@pyarts.workspace.arts_agenda
def water_psat_agenda(ws):
    ws.water_p_eq_fieldMK05()
arts.Copy(arts.water_p_eq_agenda, water_psat_agenda)
    
# %% Calculations

arts.jacobianOff()
arts.iy_unit = "PlanckBT"
arts.ppath_lmax = 10e3
arts.ppath_lraytrace = 1e3
arts.rt_integration_option = "default"
arts.rte_alonglos_v = 0.0
arts.nlteOff()

# %% Species and line absorption

arts.abs_speciesSet(species=[f"O2-Z-66-{CENTRAL_LINE_FREQ-1}-{CENTRAL_LINE_FREQ+1}"])
arts.abs_lines_per_speciesReadSpeciesSplitCatalog(basename = os.path.join(LINEPATH, ""))
arts.isotopologue_ratiosInitFromBuiltin()
arts.partition_functionsInitFromBuiltin()
arts.Wigner6Init()
                                                  
# %% Grids and planet

arts.p_grid = np.logspace(np.log10(105000), np.log10(0.1))
arts.lat_grid = np.linspace(-90, 90)
arts.lon_grid = np.linspace(-180, 180)
arts.refellipsoidEarth(model = "Sphere")
arts.z_surfaceConstantAltitude(altitude=0.0)
arts.t_surface = 293.15 + np.ones_like(arts.z_surface.value)

# %% Atmosphere

arts.AtmRawRead(basename = os.path.join(ATMPATH, ""))
arts.AtmosphereSet3D()
arts.AtmFieldsCalcExpand1D()
arts.Touch(arts.wind_u_field)
arts.Touch(arts.wind_v_field)
arts.Touch(arts.wind_w_field)
if DYNMAG:
    arts.Touch(arts.time)
    arts.MagFieldsCalcIGRF()
else:
    arts.Touch(arts.mag_u_field)
    arts.Touch(arts.mag_v_field)
    arts.Touch(arts.mag_w_field)
arts.cloudboxOff()

# %% Sensor

arts.stokes_dim = 4
arts.f_grid = np.linspace(-SPECTROMETER_HW, SPECTROMETER_HW, NF) + CENTRAL_LINE_FREQ
arts.sensor_pos = np.zeros((NR, 3))
arts.sensor_pos.value[:, 0] = 300e3
arts.sensor_pos.value[:, 1] = np.linspace(-80, 80, NR)
arts.sensor_los = np.zeros((NR, 2))
arts.sensor_los.value[:, 0] = 180
arts.sensorOff()

# %% Computechecks

arts.atmgeom_checkedCalc()
arts.lbl_checkedCalc()
arts.atmfields_checkedCalc()
arts.cloudbox_checkedCalc()
arts.sensor_checkedCalc()
arts.propmat_clearsky_agenda_checkedCalc()

# %% Compute

arts.yCalc()

# %% Save and compare data

if DO_SAVE:
    pyarts.xml.save(arts.y.value, os.path.join(testdir, "refdata.xml"), precision='g')

if CF_SAVE:
    arts.VectorCreate("yref")
    arts.ReadXML(arts.yref, os.path.join(testdir, "refdata.xml"))
    
    if SHOW_PLOTS:
        f = (arts.f_grid.value - CENTRAL_LINE_FREQ) / 1e6  # MHz
        y = arts.yref.value
        
        plt.plot(f, y[::4].reshape(NR, NF).T)
        plt.xlabel("Freq offset [MHz]")
        plt.ylabel("I [K]")
        plt.title("Downlooking at 0-longitude [REFERENCE]")
        plt.legend(arts.sensor_pos.value[:, 1], title="Latitude", loc='lower left')
        plt.show()
        
        plt.plot(f, y[1::4].reshape(NR, NF).T)
        plt.xlabel("Freq offset [MHz]")
        plt.ylabel("Q [K]")
        plt.title("Downlooking at 0-longitude [REFERENCE]")
        plt.legend(arts.sensor_pos.value[:, 1], title="Latitude", loc='lower left')
        plt.show()
        
        plt.plot(f, y[2::4].reshape(NR, NF).T)
        plt.xlabel("Freq offset [MHz]")
        plt.ylabel("U [K]")
        plt.title("Downlooking at 0-longitude [REFERENCE]")
        plt.legend(arts.sensor_pos.value[:, 1], title="Latitude", loc='lower left')
        plt.show()
        
        plt.plot(f, y[3::4].reshape(NR, NF).T)
        plt.xlabel("Freq offset [MHz]")
        plt.ylabel("V [K]")
        plt.title("Downlooking at 0-longitude [REFERENCE]")
        plt.legend(arts.sensor_pos.value[:, 1], title="Latitude", loc='lower left')
        plt.show()
    
    arts.CompareRelative(arts.yref, arts.y, 1e-5, "y reference validation failed")

#%% Plot current

if SHOW_PLOTS:
    f = (arts.f_grid.value - CENTRAL_LINE_FREQ) / 1e6  # MHz
    plt.plot(f, arts.y.value[::4].reshape(NR, NF).T)
    plt.xlabel("Freq offset [MHz]")
    plt.ylabel("I [K]")
    plt.title("Downlooking at 0-longitude")
    plt.legend(arts.sensor_pos.value[:, 1], title="Latitude", loc='lower left')
    plt.show()
    
    plt.plot(f, arts.y.value[1::4].reshape(NR, NF).T)
    plt.xlabel("Freq offset [MHz]")
    plt.ylabel("Q [K]")
    plt.title("Downlooking at 0-longitude")
    plt.legend(arts.sensor_pos.value[:, 1], title="Latitude", loc='lower left')
    plt.show()
    
    plt.plot(f, arts.y.value[2::4].reshape(NR, NF).T)
    plt.xlabel("Freq offset [MHz]")
    plt.ylabel("U [K]")
    plt.title("Downlooking at 0-longitude")
    plt.legend(arts.sensor_pos.value[:, 1], title="Latitude", loc='lower left')
    plt.show()
    
    plt.plot(f, arts.y.value[3::4].reshape(NR, NF).T)
    plt.xlabel("Freq offset [MHz]")
    plt.ylabel("V [K]")
    plt.title("Downlooking at 0-longitude")
    plt.legend(arts.sensor_pos.value[:, 1], title="Latitude", loc='lower left')
    plt.show()
