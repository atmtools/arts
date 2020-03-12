# DEFINITIONS:  -*-sh-*-
# A simple test and demonstration of *yCalcAppend*.
#
# A retrieval based on a combination of three sources is assumed:
# 1. Emission spectrum of 110.8 GHz ozone transition
# 2. A solar occultation spectrum of the same transition
# 3. Data from a 2-channel tropospheric microwave radiometer
#
# The set-up and selection of retrieval quantities are not realistic, rather
# selected to test the code in various ways.
#
# Author: Patrick Eriksson

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/planet_earth.arts")
# Agenda for scalar gas absorption calculation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
# standard surface agenda (i.e., make use of surface_rtprop_agenda)
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
# sensor-only path
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
# no refraction
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
# On-the-fly absorption
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
# ---- Basic and common stuff -----------------------------------------------
# Number of Stokes components to be computed
ws.IndexSet(ws.stokes_dim, 1)
# Frequency grid
ws.VectorNLinSpace(ws.f_grid, 201, 110000000000.0, 111000000000.0)
# Dimensionality of the atmosphere
ws.AtmosphereSet1D()
# Species and absorption
ws.abs_speciesSet(species=["O3", "H2O-PWR98", "N2-SelfContStandardType", "O2-PWR98"])
ws.ReadARTSCAT(ws.abs_lines, "ozone_line.xml")
ws.abs_lines_per_speciesCreateFromLines()
# Surface altitude
ws.MatrixSetConstant(ws.z_surface, 1, 1, 0.0)
# A pressure grid rougly matching 0 to 80 km in 500 m steps.
ws.VectorNLogSpace(ws.p_grid, 160, 101300.0, 1.0)
ws.AtmRawRead(basename="testdata/tropical")
ws.AtmFieldsCalc()
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.lbl_checkedCalc()
# ---- Part 1: An emission measurement at 110.8 GHz --------------------------
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
ws.MatrixSetConstant(ws.sensor_pos, 1, 1, 0.0)
ws.MatrixSetConstant(ws.sensor_los, 1, 1, 0.0)
ws.sensorOff()
ws.sensor_checkedCalc()
ws.StringSet(ws.iy_unit, "RJBT")
ws.ArrayOfStringSet(ws.iy_aux_vars, ["Optical depth", "Radiative background"])
ws.jacobianInit()
ws.jacobianAddAbsSpecies(
    g1=ws.p_grid, g2=ws.lat_grid, g3=ws.lon_grid, species="O3", unit="rel"
)
ws.jacobianAddPolyfit(poly_order=1)
ws.jacobianClose()
# No cloudbox
ws.cloudboxOff()
ws.cloudbox_checkedCalc()
ws.yCalc()
# ---- Part 2: Corresponding transmission (at another angle) -----------------
ws.Copy(ws.iy_transmitter_agenda, ws.iy_transmitter_agenda__UnitUnpolIntensity)
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Transmission)
ws.StringSet(ws.iy_unit, "1")
ws.MatrixSetConstant(ws.sensor_los, 1, 1, 45.0)
ws.sensorOff()
ws.sensor_checkedCalc()
ws.ArrayOfStringSet(ws.iy_aux_vars, ["Optical depth"])
ws.ArrayOfRetrievalQuantityCreate("jacobian_quantities_copy")
ws.Copy(ws.jacobian_quantities_copy, ws.jacobian_quantities)
ws.jacobianInit()
ws.jacobianAddAbsSpecies(
    g1=ws.p_grid, g2=ws.lat_grid, g3=ws.lon_grid, species="H2O-PWR98", unit="rel"
)
ws.VectorCreate("rgrid")
# Here just to test check of consistency between grids
# VectorNLogSpace( rgrid, 159, 1.013e5, 1 )
ws.Copy(ws.rgrid, ws.p_grid)
ws.jacobianAddAbsSpecies(
    g1=ws.rgrid, g2=ws.lat_grid, g3=ws.lon_grid, species="O3", unit="rel"
)
ws.jacobianAddPolyfit(poly_order=0)
ws.jacobianClose()
# No cloudbox
ws.cloudboxOff(jacobian_quantities=ws.jacobian_quantities_copy)
ws.cloudbox_checkedCalc(jacobian_quantities=ws.jacobian_quantities_copy)
ws.yCalcAppend(
    jacobian_quantities_copy=ws.jacobian_quantities_copy, append_instrument_wfs=0
)
# ---- Part 3: Data from tropospheric water vapour radiometer -----------------
ws.VectorSet(ws.f_grid, np.array([2.34e10, 3.10e10]))
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
ws.MatrixSetConstant(ws.sensor_los, 1, 1, 0.0)
ws.sensorOff()
ws.sensor_checkedCalc()
ws.StringSet(ws.iy_unit, "RJBT")
ws.ArrayOfStringSet(ws.iy_aux_vars, [])
ws.Copy(ws.jacobian_quantities_copy, ws.jacobian_quantities)
ws.jacobianInit()
ws.jacobianAddAbsSpecies(
    g1=ws.p_grid, g2=ws.lat_grid, g3=ws.lon_grid, species="H2O-PWR98", unit="rel"
)
ws.jacobianClose()
# No cloudbox
ws.cloudboxOff(jacobian_quantities=ws.jacobian_quantities_copy)
ws.cloudbox_checkedCalc(jacobian_quantities=ws.jacobian_quantities_copy)
ws.yCalcAppend(jacobian_quantities_copy=ws.jacobian_quantities_copy)
# Save reults
ws.WriteXML(ws.output_file_format, ws.y, "TestYCalcAppend.y.xml")
ws.WriteXML(ws.output_file_format, ws.y_f, "TestYCalcAppend.y_f.xml")
ws.WriteXML(ws.output_file_format, ws.y_pol, "TestYCalcAppend.y_pol.xml")
ws.WriteXML(ws.output_file_format, ws.y_pos, "TestYCalcAppend.y_pos.xml")
ws.WriteXML(ws.output_file_format, ws.y_los, "TestYCalcAppend.y_los.xml")
ws.WriteXML(ws.output_file_format, ws.y_aux, "TestYCalcAppend.y_aux.xml")
#
ws.WriteXML(ws.output_file_format, ws.jacobian, "TestYCalcAppend.jacobian.xml")
ws.WriteXML(
    ws.output_file_format,
    ws.jacobian_quantities,
    "TestYCalcAppend.jacobian_quantities.xml",
)
