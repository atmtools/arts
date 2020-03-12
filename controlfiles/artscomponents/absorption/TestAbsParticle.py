# DEFINITIONS:  -*-sh-*-
#
# Demonstration of a calculation with absorption-only particles (handling them
# as abs_species).
#
# Copy of the TestDOIT case, just replacing DOIT scattering parts by
# particles-as-abs_species setup.
#
# Author: Jana Mendrok, Claudia Emde
#

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
# Agenda for scalar gas absorption calculation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
# (standard) emission calculation
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
# cosmic background radiation
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
# standard surface agenda (i.e., make use of surface_rtprop_agenda)
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
# sensor-only path
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
# no refraction
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
# No jacobian calculations
ws.jacobianOff()
# Frequency grid
# --------------
# Note: The frequencies must be contained in the gas absorption lookup table.
ws.VectorSet(ws.f_grid, np.array([2.295e11, 2.305e11]))
# Number of Stokes components to be computed
# -------------------------------------------
ws.IndexSet(ws.stokes_dim, 4)
# Definition of the atmosphere
# ----------------------------
# Dimensionality of the atmosphere
ws.AtmosphereSet1D()
# Pressure grid
ws.ReadXML(ws.p_grid, "testdata/testdoit_p_grid.xml")
# Definition of species
ws.abs_speciesSet(species=["H2O-PWR98", "O2-PWR93", "N2-SelfContStandardType"])
# Atmospheric profiles
ws.AtmRawRead(basename="testdata/tropical")
# setting cloudbox_on (to off) here already. avoids writing & re-reading
# scat_data later on (since cloudboxOff wipes out all cloud-related data).
ws.cloudboxOff()
# Specification of cloud
# -----------------------
ws.ScatSpeciesInit()
# Adding scattering elements.
# Here actually both added elements are indentical. however, for testing and for
# demonstration purposed, having 2 elements is better.
ws.ScatElementsToabs_speciesAdd(
    scat_data_files=[
        "testdata/scatData/azi-random_f229-231T214-225r100NP-1ar1_5ice.xml",
        "testdata/scatData/azi-random_f229-231T214-225r100NP-1ar1_5ice.xml",
    ],
    pnd_field_files=[
        "testdata/testdoit_pnd_field_1D.xml",
        "testdata/testdoit_pnd_field_1D.xml",
    ],
)
ws.scat_dataCalc()
ws.scat_dataCheck()
# Gas absorption from lookup table
# ---------------------------------
# for how to create lookup tables, see ../absorption/TestAbs.arts
ws.ReadXML(ws.abs_lookup, "testdata/testdoit_gas_abs_lookup.xml")
ws.abs_lookupAdapt()
# absorption from LUT
# Copy( propmat_clearsky_agenda, propmat_clearsky_agenda__LookUpTable )
@arts_agenda
def propmat_clearsky_agenda(ws):
    ws.Ignore(ws.rtp_mag)
    ws.Ignore(ws.rtp_los)
    ws.Ignore(ws.rtp_nlte)
    ws.propmat_clearskyInit()
    ws.propmat_clearskyAddFromLookup()
    ws.propmat_clearskyAddParticles()


ws.propmat_clearsky_agenda = propmat_clearsky_agenda

ws.AtmFieldsCalc()
# Definition of Earth surface
# ----------------------------
ws.MatrixSetConstant(ws.z_surface, 1, 1, 500.0)
# Properties of surface:
# - surface reflectivity from Liebe model
# - surface skin temperature interpolated from atmospheric t_field
#
ws.VectorCreate("n_t_grid")


@arts_agenda
def surface_rtprop_agenda(ws):
    ws.specular_losCalc()
    ws.InterpAtmFieldToPosition(out=ws.surface_skin_t, field=ws.t_field)
    ws.VectorSetConstant(ws.n_t_grid, 1, ws.surface_skin_t)
    ws.complex_refr_indexWaterLiebe93(
        complex_refr_index=ws.surface_complex_refr_index,
        data_f_grid=ws.f_grid,
        data_T_grid=ws.n_t_grid,
    )
    ws.surfaceFlatRefractiveIndex()


ws.surface_rtprop_agenda = surface_rtprop_agenda

# Definition of sensor position and LOS
# --------------------------------------
# This file holds the viewing angles of the sensor:
ws.IndexSet(ws.nelem, 1)
ws.VectorCreate("vector_1")
ws.VectorSetConstant(ws.vector_1, ws.nelem, 99.7841941981)
# IndexSet( nelem, 19 )
# VectorNLinSpace( vector_1, 0, 180 )
# Sensor altitude from earth surface
ws.nelemGet(ws.nelem, ws.vector_1)
ws.VectorCreate("vector_2")
ws.VectorSetConstant(ws.vector_2, ws.nelem, 95000.1)
ws.Matrix1ColFromVector(ws.sensor_pos, ws.vector_2)
ws.Matrix1ColFromVector(ws.sensor_los, ws.vector_1)
# SensorOff means that the result of the calculation are the radiances,
# which are not modified by sensor properties
ws.sensorOff()
# ==================start==========================
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.scat_data_checkedCalc()
ws.sensor_checkedCalc()
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
# WriteXML( in=abs_species )
ws.StringSet(ws.iy_unit, "RJBT")
ws.yCalc()
ws.WriteXML(ws.output_file_format, ws.y, "", 0)
# ==================stop==========================
# ==================check==========================
ws.VectorCreate("yREFERENCE")
ws.ReadXML(ws.yREFERENCE, "yREFERENCE_AbsParticle.xml")
ws.Compare(ws.y, ws.yREFERENCE, 1e-06)
# End of Main
