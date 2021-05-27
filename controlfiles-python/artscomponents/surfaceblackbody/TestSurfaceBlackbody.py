import numpy as np
from pyarts.workspace import Workspace, arts_agenda


def setup_testcase(ws):
    ws.execute_controlfile("general/general.arts")
    ws.execute_controlfile("general/continua.arts")
    ws.execute_controlfile("general/agendas.arts")
    ws.execute_controlfile("general/planet_earth.arts")
    # (standard) emission calculation
    ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
    # cosmic background radiation
    ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)

    @arts_agenda
    def iy_surface_agenda_PY(ws):
        ws.SurfaceBlackbody()
        ws.iySurfaceRtpropCalc()
    # upwelling intensity (and jacobian) from the surface for given point and direction
    ws.Copy(ws.iy_surface_agenda, iy_surface_agenda_PY)
    # called by iySurfaceRtpropCalc within iy_surface_agenda
    ws.Copy(ws.surface_rtprop_agenda, ws.surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_surface)
    # clearsky agenda
    ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
    # sensor-only path
    ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
    # no refraction
    ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
    ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)

    ws.stokes_dim = 1
    ws.atmosphere_dim = 1  # 1D VAR
    ws.jacobian_quantities = []
    ws.iy_unit = "PlanckBT"
    ws.cloudboxOff()
    ws.surface_scalar_reflectivity = np.array([0.00])  # nominal albedo for surface

    # Sensor settings
    ws.sensor_pos = np.array([[850e3]])  # 850km
    ws.sensor_los = np.array([[180.]])  # nadir viewing
    ws.f_backend = np.array([183e9])
    ws.VectorNLinSpace(ws.f_grid, 10, 173e9, 193e9)
    ws.backend_channel_responseGaussian(fwhm=np.array([5e9]))
    ws.FlagOn(ws.sensor_norm)
    ws.AntennaOff()
    ws.sensor_responseInit()
    ws.sensor_responseBackend()

    ws.abs_speciesSet(species=[
        "H2O, H2O-SelfContCKDMT252, H2O-ForeignContCKDMT252",
    ])
    ws.abs_lines_per_speciesReadSpeciesSplitCatalog(
        basename="spectroscopy/Artscat/"
    )
    ws.abs_lines_per_speciesSetCutoff(option="ByLine", value=750e9)
    ws.abs_lines_per_speciesCompact()
    ws.abs_xsec_agenda_checkedCalc()
    ws.lbl_checkedCalc()

    # Load atmospheric data to be forward simulated
    ws.ReadXML(ws.batch_atm_fields_compact, 'testdata/garand_profiles.xml.gz')
    ws.propmat_clearsky_agenda_checkedCalc()

    ws.VectorCreate("t_surface_vector")  # helper variable for ybatch_calc_agenda
    ws.NumericCreate("t_surface_numeric")  # helper variable for ybatch_calc_agenda
    @arts_agenda
    def ybatch_calc_agenda(ws):
        ws.Extract(ws.atm_fields_compact,
                   ws.batch_atm_fields_compact,
                   ws.ybatch_index)
        ws.AtmFieldsAndParticleBulkPropFieldFromCompact()
        ws.Extract(ws.z_surface, ws.z_field, 0)
        ws.Extract(ws.t_surface, ws.t_field, 0)
        ws.Copy(ws.surface_props_names, ["Skin temperature"])
        ws.VectorExtractFromMatrix(ws.t_surface_vector, ws.t_surface, 0, "row")
        ws.Extract(ws.t_surface_numeric, ws.t_surface_vector, 0)
        ws.Tensor3SetConstant(ws.surface_props_data, 1, 1, 1, ws.t_surface_numeric)
        ws.jacobianInit()
        ws.jacobianAddSurfaceQuantity(
            g1=ws.lat_grid, g2=ws.lon_grid, quantity="Skin temperature")
        ws.jacobianClose()
        ws.atmfields_checkedCalc()
        ws.atmgeom_checkedCalc()
        ws.cloudbox_checkedCalc()
        ws.sensor_checkedCalc()
        ws.yCalc()

    ws.Copy(ws.ybatch_calc_agenda, ybatch_calc_agenda)
    ws.IndexSet(ws.ybatch_start, 0)
    ws.IndexSet(ws.ybatch_n, 1)  # Amount of atmospheres
    ws.ybatchCalc()  # conduct the forward simulation
    return ws


# test for ybatch
def test_ybatch(ws):
    ybatch_ref = np.array([256.96294508])
    assert np.allclose(ws.ybatch.value[0], ybatch_ref)


# test for ybatch_jacobians
def test_ybatch_jacobians(ws):
    ybatch_jacobians_ref = np.array([[3.38698132e-6]])
    assert np.allclose(ws.ybatch_jacobians.value[0], ybatch_jacobians_ref, atol=1e-12)


if __name__ == '__main__':
    ws = Workspace(verbosity=2)
    ws = setup_testcase(ws)
    test_ybatch(ws)
    test_ybatch_jacobians(ws)
    print('SurfaceBlackbody tests passed.')
