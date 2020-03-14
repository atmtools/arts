import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")


@arts_agenda
def propmat_clearsky_agenda(ws):
    ws.Ignore(ws.rtp_mag)
    ws.Ignore(ws.rtp_los)
    ws.propmat_clearskyInit()
    ws.propmat_clearskyAddOnTheFly()


ws.propmat_clearsky_agenda = propmat_clearsky_agenda


@arts_agenda
def abs_xsec_agenda(ws):
    ws.abs_xsec_per_speciesInit()
    ws.abs_xsec_per_speciesAddLines()


ws.abs_xsec_agenda = abs_xsec_agenda

ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
ws.Copy(ws.surface_rtprop_agenda, ws.surface_rtprop_agenda__Blackbody_SurfTFromt_field)
ws.IndexSet(ws.stokes_dim, 1)
ws.VectorNLogSpace(ws.p_grid, 41, 100000.0, 1.0)
ws.VectorNLinSpace(ws.lat_grid, 2, -60.0, 60.0)
ws.VectorNLinSpace(ws.lon_grid, 2, -30.0, 30.0)
ws.MatrixSet(
    ws.sensor_los,
    np.array([[112.68855143, 0.0], [112.25892819, 0.0], [111.82133233, 0.0]]),
)
ws.MatrixSet(
    ws.sensor_pos,
    np.array([[600000.0, 0.0, 0.0], [600000.0, 0.0, 0.0], [600000.0, 0.0, 0.0]]),
)
ws.refellipsoidEarth(ws.refellipsoid, "Sphere")
ws.AtmosphereSet3D()
ws.abs_speciesSet(species=["CO2-626"])
ws.VectorCreate("tmp_kayser")
ws.VectorLinSpace(ws.tmp_kayser, 600.0, 650.0, 1.0)
ws.FrequencyFromCGSKayserWavenumber(ws.f_grid, ws.tmp_kayser)
# WriteXML("ascii",tmp_kayser,"cm.xml")
# WriteXML("ascii",f_grid,"f.xml")
ws.ReadArrayOfARTSCAT(
    abs_lines=ws.abs_lines,
    filename="../../testdata/NLTE_CO2_testlines.xml",
    globalquantumnumbers="v1 v2 v3 r l2",
    localquantumnumbers="J",
)
ws.abs_linesSetNormalization(option="VVH")
ws.abs_linesSetCutoff(option="ByLine", value=750000000000.0)
ws.abs_lines_per_speciesCreateFromLines()
ws.isotopologue_ratiosInitFromBuiltin()
ws.sensorOff()
ws.jacobianOff()
ws.cloudboxOff()
ws.StringSet(ws.iy_unit, "W/(m^2 m-1 sr)")
# Note that this data is not at all reliable and does not represent a real world scenario...
ws.AtmWithNLTERawRead(basename="testdata/tropical", expect_vibrational_energies=1)
ws.AtmFieldsCalcExpand1D()
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.VectorCreate("ev")
ws.ReadXML(ws.ev, "testdata/tropical.ev.xml")
ws.nlteSetByQuantumIdentifiers()
# NumericSet(ppath_lmax,5e3)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
ws.lbl_checkedCalc()
ws.yCalc()
ws.VectorCreate("y_nlte_ref")
ws.WriteXML("ascii", ws.y, "TestNLTE_NLTE.xml")
ws.ReadXML(ws.y_nlte_ref, "TestNLTE_NLTE_REFERENCE.xml")
ws.CompareRelative(ws.y, ws.y_nlte_ref, 1e-05, "NLTE old method")
ws.EnergyLevelMapCreate("tmp_elm")
ws.ArrayOfQuantumIdentifierCreate("tmp_aoqi")
ws.Copy(ws.tmp_elm, ws.nlte_field)
ws.Copy(ws.tmp_aoqi, ws.nlte_level_identifiers)
ws.nlteOff()
ws.abs_lines_per_speciesSetPopulation(option="LTE")
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
ws.yCalc()
ws.VectorCreate("y_lte_ref")
ws.WriteXML("ascii", ws.y, "TestNLTE_LTE.xml")
ws.ReadXML(ws.y_lte_ref, "TestNLTE_LTE_REFERENCE.xml")
ws.CompareRelative(ws.y, ws.y_lte_ref, 1e-05, "LTE old method")
