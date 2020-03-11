# DEFINITIONS:  -*-sh-*-
# This control file performs an ARTS-MC radiative transfer simulation
# For a single line of sight and a Gaussian antenna response

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
ws.jacobianOff()
# Agenda for scalar gas absorption calculation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
# cosmic background radiation
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
# no refraction
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
# blackbody surface with skin temperature interpolated from t_surface field
ws.Copy(ws.surface_rtprop_agenda, ws.surface_rtprop_agenda__Blackbody_SurfTFromt_field)
#### LOAD DATA: these files were created with MCDataPrepare.arts ######
ws.ReadXML(ws.f_grid, "TestMonteCarloDataPrepare.f_grid.xml")
ws.IndexSet(ws.f_index, 0)
ws.ReadXML(ws.p_grid, "p_grid.xml")
ws.AtmosphereSet3D()
ws.ReadXML(ws.lat_grid, "lat_grid.xml")
ws.ReadXML(ws.lon_grid, "lon_grid.xml")
ws.ReadXML(ws.t_field, "TestMonteCarloDataPrepare.t_field.xml")
ws.ReadXML(ws.z_field, "TestMonteCarloDataPrepare.z_field.xml")
ws.ReadXML(ws.vmr_field, "TestMonteCarloDataPrepare.vmr_field.xml")
ws.ReadXML(ws.z_surface, "TestMonteCarloDataPrepare.z_surface.xml")
ws.ReadXML(ws.abs_lookup, "TestMonteCarloDataPrepare.abs_lookup.xml")
ws.abs_speciesSet(species=["O2-PWR93", "N2-SelfContStandardType", "H2O-PWR98"])
ws.abs_lookupAdapt()
ws.FlagOn(ws.cloudbox_on)
ws.ReadXML(ws.cloudbox_limits, "TestMonteCarloDataPrepare.cloudbox_limits.xml")
ws.ReadXML(ws.pnd_field, "TestMonteCarloDataPrepare.pnd_field.xml")
ws.ReadXML(ws.scat_data, "TestMonteCarloDataPrepare.scat_data.xml")
ws.scat_data_checkedCalc()
#### Define Agendas #################################################
# absorption from LUT
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__LookUpTable)
#### Define viewing position and line of sight #########################
ws.rte_losSet(ws.rte_los, ws.atmosphere_dim, 99.7841941981, 180.0)
ws.rte_posSet(ws.rte_pos, ws.atmosphere_dim, 95000.1, 7.61968838781, 0.0)
ws.Matrix1RowFromVector(ws.sensor_pos, ws.rte_pos)
ws.Matrix1RowFromVector(ws.sensor_los, ws.rte_los)
#### Set some Monte Carlo parameters ###################################
ws.IndexSet(ws.stokes_dim, 4)
ws.StringSet(ws.iy_unit, "RJBT")
ws.NumericSet(ws.ppath_lmax, 3000.0)
ws.MCSetSeedFromTime()
ws.Print(ws.mc_seed, 1)
ws.mc_antennaSetGaussianByFWHM(ws.mc_antenna, 0.1137, 0.239)
#### Check atmosphere ##################################################
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
#### Perform Monte Carlo RT Calculation #################################
ws.NumericSet(ws.mc_std_err, -1.0)
ws.IndexSet(ws.mc_max_time, 20)
ws.IndexSet(ws.mc_max_iter, -1)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.MCGeneral()
#### Save calculated Stokes vector and std. err.#########################
ws.output_file_formatSetAscii()
ws.WriteXML(ws.output_file_format, ws.y)
ws.WriteXML(ws.output_file_format, ws.mc_error)
#### Print number of photons and radiance units ########################
ws.Print(ws.mc_iteration_count, 1)
ws.Print(ws.iy_unit, 1)
#### Tests ########################
# Radiance test
ws.NumericCreate("mc_error_0")
ws.Extract(ws.mc_error_0, ws.mc_error, 0)
ws.NumericScale(ws.mc_error_0, ws.mc_error_0, 4.0)
ws.VectorCreate("y_0")
ws.Select(ws.y_0, ws.y, [0])
ws.VectorCreate("y_ref")
ws.VectorSetConstant(ws.y_ref, 1, 198.6)
ws.Compare(ws.y_0, ws.y_ref, ws.mc_error_0, "Total radiance should be close to 198.6 K")
# Polarization test
ws.NumericCreate("mc_error_1")
ws.Extract(ws.mc_error_1, ws.mc_error, 1)
ws.NumericScale(ws.mc_error_1, ws.mc_error_1, 4.0)
ws.VectorCreate("y_1")
ws.Select(ws.y_1, ws.y, [1])
ws.VectorSetConstant(ws.y_ref, 1, 7.6)
ws.Compare(
    ws.y_1, ws.y_ref, ws.mc_error_1, "Polarization difference should be close to 7.6 K"
)
