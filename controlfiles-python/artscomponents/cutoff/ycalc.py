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

ATMPATH = "../zeeman/atm"
LINEPATH = "../zeeman/line"

SHOW_PLOTS = False

DO_SAVE = False

CF_SAVE = True

# %% Optional imports

if SHOW_PLOTS:
    import matplotlib.pyplot as plt

# %% Silly ARTS-isms
    
arts.abs_f_interp_order = 100_000_000
arts.Touch(arts.iy_aux_vars)
arts.Touch(arts.surface_props_data)
arts.Touch(arts.surface_props_names)
arts.Touch(arts.transmitter_pos)

# %% Agendas

@pyarts.workspace.arts_agenda(ws=arts)
def propmat_clearsky_agenda(arts):
    arts.propmat_clearskyInit()
    arts.propmat_clearskyAddLines()
    arts.Ignore(arts.rtp_mag)
    arts.Ignore(arts.rtp_los)
arts.Copy(arts.propmat_clearsky_agenda, propmat_clearsky_agenda)

@pyarts.workspace.arts_agenda(ws=arts)
def ppath_agenda_step_by_step(arts):
    arts.Ignore(arts.rte_pos2)
    arts.ppathStepByStep()
arts.Copy(arts.ppath_agenda, ppath_agenda_step_by_step)

@pyarts.workspace.arts_agenda(ws=arts)
def iy_main_agenda_emission(arts):
    arts.ppathCalc()
    arts.iyEmissionStandard()
arts.Copy(arts.iy_main_agenda, iy_main_agenda_emission)

@pyarts.workspace.arts_agenda(ws=arts)
def surface_rtprop_agenda(arts):
    arts.InterpSurfaceFieldToPosition(out=arts.surface_skin_t, field=arts.t_surface)
    arts.surfaceBlackbody()
arts.Copy(arts.surface_rtprop_agenda, surface_rtprop_agenda)
 
@pyarts.workspace.arts_agenda(ws=arts)
def ppath_step_agenda_geometric(arts):
    arts.Ignore(arts.t_field)
    arts.Ignore(arts.vmr_field)
    arts.Ignore(arts.f_grid)
    arts.Ignore(arts.ppath_lraytrace)
    arts.ppath_stepGeometric()
arts.Copy(arts.ppath_step_agenda, ppath_step_agenda_geometric)

@pyarts.workspace.arts_agenda(ws=arts)
def iy_space_agenda_cosmic_background(arts):
    arts.Ignore(arts.rtp_pos)
    arts.Ignore(arts.rtp_los)
    arts.MatrixCBR(arts.iy, arts.stokes_dim, arts.f_grid)
arts.Copy(arts.iy_space_agenda, iy_space_agenda_cosmic_background)

@pyarts.workspace.arts_agenda(ws=arts)
def geo_pos_agenda(arts):
    arts.Ignore(arts.ppath)
    arts.VectorSet(arts.geo_pos, np.array([]))
arts.Copy(arts.geo_pos_agenda, geo_pos_agenda)

@pyarts.workspace.arts_agenda(ws=arts)
def iy_surface_agenda(arts):
    arts.SurfaceDummy()
    arts.iySurfaceRtpropAgenda()
arts.Copy(arts.iy_surface_agenda, iy_surface_agenda)

@pyarts.workspace.arts_agenda(ws=arts)
def water_psat_agenda(arts):
    arts.water_p_eq_fieldMK05()
arts.Copy(arts.water_p_eq_agenda, water_psat_agenda)
    
# %% Calculations

arts.jacobianOff()
arts.iy_unit = "1"
arts.ppath_lmax = 10e3
arts.ppath_lraytrace = 1e3
arts.rt_integration_option = "default"
arts.rte_alonglos_v = 0.0
arts.nlteOff()

# %% Species and line absorption

arts.abs_speciesSet(species=["O2-66"])
arts.abs_lines_per_speciesReadSpeciesSplitCatalog(basename = os.path.join(LINEPATH, ""))
arts.isotopologue_ratiosInitFromBuiltin()

# %% Introduce a weird shift so that the cutoff is tested properly

x = arts.abs_lines_per_species.value[0][0].lines[0].lineshape[0].D0 = \
    pyarts.classes.LineShapeModelParameters("T0", 50000)
                                                  
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
arts.Touch(arts.mag_u_field)
arts.Touch(arts.mag_v_field)
arts.Touch(arts.mag_w_field)
arts.cloudboxOff()

# %% Sensor

arts.stokes_dim = 1
arts.f_grid = np.linspace(-5e9, 5e9, 101) + 118750348044.712
arts.sensor_pos = np.array([[300e3, 0, 0]])
arts.sensor_los = np.array([[180, 0]])
arts.sensorOff()

# %% Computechecks

arts.atmgeom_checkedCalc()
arts.lbl_checkedCalc()
arts.atmfields_checkedCalc()
arts.cloudbox_checkedCalc()
arts.sensor_checkedCalc()
arts.propmat_clearsky_agenda_checkedCalc()

# %% Compute

cutoffs = np.logspace(7, 13, 13-7+1)
y = []
for cutoff in cutoffs:
    arts.abs_lines_per_speciesSetCutoff(option="ByLine", value=cutoff)
    arts.yCalc()
    y.append(arts.y.value * 1.0)
y = np.array(y)

# %% Save and compare data

if DO_SAVE:
    pyarts.xml.save(y, os.path.join(testdir, "refdata.xml"), precision='g')

if CF_SAVE:
    arts.MatrixCreate("yref")
    arts.ReadXML(arts.yref, "refdata.xml")
    
    arts.CompareRelative(arts.yref, y, 1e-5, "y reference validation failed")

if SHOW_PLOTS:
    plt.figure(figsize=[8, 4])
    plt.title("Weirdly shifted absorption line")
    plt.semilogy(arts.f_grid.value / 1e9, y.T)
    plt.ylabel("Radiance [W/(m$^2$ Hz sr)]")
    plt.xlabel("Frequency [GHz]")
    plt.legend([f"Cutoff at {c/1e9} GHz" for c in cutoffs])
    plt.show()
    