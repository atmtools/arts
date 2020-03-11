# DEFINITIONS:  -*-sh-*-
#
# This control file performs an ARTS-MC radiative transfer simulation
# through yCalc

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
# Agenda for scalar gas absorption calculation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
# cosmic background radiation
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
# no refraction
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
# Blackbody surface with skin temperature interpolated from t_surface field
ws.Copy(ws.surface_rtprop_agenda, ws.surface_rtprop_agenda__Blackbody_SurfTFromt_field)
ws.IndexSet(ws.stokes_dim, 1)
ws.AtmosphereSet3D()
ws.jacobianOff()
#### LOAD DATA: these files were created with MCDataPrepare.arts ######
ws.ReadXML(ws.f_grid, "TestMonteCarloDataPrepare.f_grid.xml")
ws.ReadXML(ws.p_grid, "p_grid.xml")
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
ws.scat_dataCheck()
ws.scat_data_checkedCalc()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
#### Define Agendas #################################################
# absorption from LUT
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__LookUpTable)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
#### Define viewing position and line of sight #########################
ws.rte_losSet(ws.rte_los, ws.atmosphere_dim, 99.7841941981, 180.0)
ws.rte_posSet(ws.rte_pos, ws.atmosphere_dim, 95000.1, 7.61968838781, 0.0)
ws.Matrix1RowFromVector(ws.sensor_pos, ws.rte_pos)
ws.Matrix1RowFromVector(ws.sensor_los, ws.rte_los)
#### Include an antenna pattern #########################################
ws.IndexSet(ws.stokes_dim, 4)
ws.IndexSet(ws.antenna_dim, 1)
#
ws.VectorSet(ws.za_grid, array([-0.2, -0.02, 0.0, 0.02, 0.2]))
ws.Matrix1ColFromVector(ws.mblock_dlos_grid, ws.za_grid)
#
ws.MatrixSetConstant(ws.antenna_dlos, 1, 1, 0.0)
# An antenna pattern from Odin-SMR is used here
ws.ReadXML(ws.antenna_response, "antenna.SM_AC2ab.875ms.xml")
ws.IndexSet(ws.sensor_norm, 1)
ws.sensor_responseInit()
ws.sensor_responseAntenna()
ws.sensor_checkedCalc()
#### Perform Monte Carlo RT Calculation #################################
ws.StringSet(ws.iy_unit, "RJBT")
ws.NumericSet(ws.ppath_lmax, 3000.0)
ws.NumericSet(ws.mc_std_err, -1.0)
ws.IndexSet(ws.mc_max_time, -1)
ws.IndexSet(ws.mc_max_iter, 100)
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__ScattMC)
ws.ArrayOfStringSet(ws.iy_aux_vars, ["Error (uncorrelated)"])
ws.yCalc()
#### Save results #########################
# output_file_formatSetAscii
# WriteXML( output_file_format, y )
# WriteXML( output_file_format, y_aux )
#### Tests ########################
# Radiance test
ws.VectorCreate("vtmp")
ws.Extract(ws.vtmp, ws.y_aux, 0)
ws.NumericCreate("y_error_0")
ws.Extract(ws.y_error_0, ws.vtmp, 0)
ws.NumericScale(ws.y_error_0, ws.y_error_0, 4.0)
ws.VectorCreate("y_0")
ws.Select(ws.y_0, ws.y, [0])
ws.VectorCreate("y_ref")
ws.VectorSetConstant(ws.y_ref, 1, 198.7)
ws.Print(ws.y_0, 0)
ws.Print(ws.y_ref, 0)
ws.Print(ws.y_error_0, 0)
ws.Compare(ws.y_0, ws.y_ref, ws.y_error_0, "Total radiance should be close to 198.7")
