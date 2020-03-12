# DEFINITIONS:  -*-sh-*-
#
# A simple test of yApplySensorPol. The test just checks that teh code runs,
# there is no comparison to some reference values.
#
# 2015-10-11, Patrick Eriksson

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)


@arts_agenda
def geo_pos_agenda(ws):
    ws.geo_posEndOfPpath()


ws.geo_pos_agenda = geo_pos_agenda

#
# Create an empty (having zero absorption) 3D atmosphere
#
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
ws.abs_speciesSet(species=["H2O", "N2"])
ws.abs_lines_per_speciesSetEmpty()
ws.VectorSet(ws.p_grid, np.array([1.013e05, 1.000e00]))
ws.VectorSet(ws.lat_grid, np.array([-90.0, 90.0]))
ws.VectorSet(ws.lon_grid, np.array([-180.0, 180.0]))
ws.AtmosphereSet3D()
ws.MatrixSetConstant(ws.z_surface, 2, 2, 0.0)
ws.Tensor3SetConstant(ws.t_field, 2, 2, 2, 300.0)
ws.Tensor4SetConstant(ws.vmr_field, 2, 2, 2, 2, 0.0)
ws.Tensor3SetConstant(ws.z_field, 2, 2, 2, 0.0)
ws.atmfields_checkedCalc()
ws.z_fieldFromHSE(p_hse=101300.0, z_hse_accuracy=1000.0)
#
# A flat water surface, at 300K
#
ws.VectorCreate("data_t_grid")
ws.VectorSet(ws.data_t_grid, np.array([290.0, 300.0, 310.0]))
ws.VectorCreate("data_f_grid")
ws.VectorLinSpace(ws.data_f_grid, 10000000000.0, 100000000000.0, 5000000000.0)
ws.complex_refr_indexWaterLiebe93(
    ws.surface_complex_refr_index, ws.data_f_grid, ws.data_t_grid
)
#
@arts_agenda
def surface_rtprop_agenda(ws):
    ws.NumericSet(ws.surface_skin_t, 300.0)
    ws.specular_losCalc()
    ws.surfaceFlatRefractiveIndex()


ws.surface_rtprop_agenda = surface_rtprop_agenda

# Number of Stokes components to be computed
#
ws.IndexSet(ws.stokes_dim, 3)
# Frequency grid
#
ws.VectorSet(ws.f_grid, np.array([1.0e10, 2.0e10, 3.0e10]))
# Sensor pos, los and pol
#
ws.MatrixSet(
    ws.sensor_pos, np.array([[6.0e05, 2.3e01, 4.0e00], [6.0e05, 7.8e01, 7.7e01]])
)
ws.MatrixSet(ws.sensor_los, np.array([[145.0, 30.0], [145.0, -26.0]]))
ws.MatrixSet(ws.sensor_pol, np.array([[0.0, 90.0, -45.0], [0.0, 90.0, -45.0]]))
# No "standard" sensor responses, but use RJ-Tb
#
ws.sensorOff()
ws.StringSet(ws.iy_unit, "RJBT")
# Add some Jacobians, to test also thos part
ws.jacobianInit()
ws.nrowsGet(ws.nrows, ws.sensor_pos)
ws.VectorNLinSpace(ws.sensor_time, ws.nrows, 0.0, 1.0)
ws.jacobianAddPointingZa(
    ws.jacobian_quantities,
    ws.jacobian_agenda,
    ws.sensor_pos,
    ws.sensor_time,
    0,
    "recalc",
    0.001,
)
ws.jacobianClose()
# No scattering
#
ws.cloudboxOff()
# Perform remaining checks
#
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
ws.lbl_checkedCalc()
# Calculate Stokes vector values
ws.ArrayOfStringSet(ws.iy_aux_vars, ["Optical depth", "Radiative background"])
ws.yCalc()
# Print( y )
# Print( jacobian )
# Extrat Tb for selected polarisation angles
ws.yApplySensorPol()
# Print( y )
# Print( jacobian )
