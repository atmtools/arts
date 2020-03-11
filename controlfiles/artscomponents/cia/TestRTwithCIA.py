# DEFINITIONS:  -*-sh-*-
#
# Demonstration and test the use of CIA contiua in RT calculations,
# together with other absorption tags.
#
# Based on Patrick Eriksson's TestClearSky2.arts
#
# Author: Stefan Buehler
# Date:   2013-03-08

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
# (standard) emission calculation
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
# cosmic background radiation
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
# standard surface agenda (i.e., make use of surface_rtprop_agenda)
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
# on-the-fly absorption
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
# Include CIA in spectroscopy
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__withCIA)
# sensor-only path
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
# no refraction
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
# Number of Stokes components to be computed
#
ws.IndexSet(ws.stokes_dim, 1)
# No jacobian calculation
#
ws.jacobianOff()
# Clearsky = No scattering
#
ws.cloudboxOff()
# Definition of species
# ---
# ATTENTION: You also have to include CIA in abs_xsec_agenda,
# otherwise you will not get any absorption for it!
ws.abs_speciesSet(
    species=[
        "H2O-SelfContCKDMT100, H2O-ForeignContCKDMT100, H2O",
        "O2",
        "N2-CIA-N2-0, N2-CIA-N2-1",
    ]
)
# Read CIA data
# ---
ws.ReadXML(ws.abs_cia_data, "spectroscopy/cia/hitran2011/hitran_cia2012_adapted.xml.gz")
# Line file and a matching coarse frequency grid
# ---
ws.ReadSplitARTSCAT(
    basename="spectroscopy/Perrin/", fmin=1000000000.0, fmax=3000000000000.0
)
ws.abs_linesSetCutoff(option="ByLine", value=750000000000.0)
ws.abs_linesSetNormalization(option="VVH")
ws.VectorNLinSpace(ws.f_grid, 100, 1000000000.0, 3000000000000.0)
# A pressure grid rougly matching 0 to 80 km, in steps of 2 km.
# ---
ws.VectorNLogSpace(ws.p_grid, 41, 100000.0, 1.0)
# Sort the line file according to species
# ---
ws.abs_lines_per_speciesCreateFromLines()
# Atmospheric scenario
# ---
ws.AtmRawRead(basename="testdata/tropical")
# Weakly reflecting surface
# ---
ws.VectorSetConstant(ws.surface_scalar_reflectivity, 1, 0.8)
ws.Copy(
    ws.surface_rtprop_agenda,
    ws.surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_surface,
)
# Tb unit
# ---
ws.StringSet(ws.iy_unit, "PlanckBT")
# Atmosphere and surface
# ---
ws.AtmosphereSet1D()
ws.AtmFieldsCalc()
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.Extract(ws.t_surface, ws.t_field, 0)
# Definition of sensor position and LOS
# ---
ws.VectorSet(ws.rte_pos, array([600000.0]))
ws.VectorSet(ws.rte_los, array([180.0]))
ws.VectorSet(ws.rte_pos2, array([], dtype=float64))
# A dummy value
# Perform RT calculations
# ---
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.lbl_checkedCalc()
ws.iyCalc()
# WriteXML( "ascii", iy )
# WriteXML( "ascii", f_grid )
# WriteXML( "ascii", abs_species )
# Compare results to reference calculation:
ws.MatrixCreate("iy_reference")
ws.ReadXML(ws.iy_reference, "TestRTwithCIA.iy_reference.xml")
ws.Compare(
    ws.iy,
    ws.iy_reference,
    1e-08,
    "Brightness temperature with CIA N2 continuum differs from reference calculation",
)
