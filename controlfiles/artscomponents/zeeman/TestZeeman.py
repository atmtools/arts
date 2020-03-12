import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.Wigner6Init(ws.wigner_initialized, 40000, 100)
# Set lines and species
ws.abs_speciesSet(species=["O2-Z-66"])
ws.ReadARTSCAT(
    abs_lines=ws.abs_lines,
    filename="testdata/zeeman-lines.xml",
    localquantumnumbers="J",
)
ws.abs_lines_per_speciesCreateFromLines()
# Initialize standard inputs
ws.isotopologue_ratiosInitFromBuiltin()
ws.partition_functionsInitFromBuiltin()
ws.VectorNLinSpace(ws.f_grid, 501, 99990000000.0, 100010000000.0)
ws.Touch(ws.rtp_nlte)
ws.IndexSet(ws.propmat_clearsky_agenda_checked, 1)
ws.AtmosphereSet3D()
ws.IndexSet(ws.stokes_dim, 4)
ws.IndexSet(ws.nlte_do, 0)
ws.IndexSet(ws.abs_f_interp_order, 3)
# Set atmospheric variables
ws.NumericSet(ws.rtp_pressure, 10.0)
ws.NumericSet(ws.rtp_temperature, 215.0)
ws.VectorSet(ws.rtp_vmr, np.array([0.21]))
ws.VectorSet(ws.rtp_mag, np.array([2.5e-05, 6.0e-05, 1.0e-05]))
ws.VectorSet(ws.rtp_los, np.array([60.0, 50.0]))
# Test variables
ws.ArrayOfPropagationMatrixCreate("test")
## Silly parameters that have to be set by agendas and ARTS in general but are completely useless for these calculations
ws.VectorSet(ws.p_grid, np.array([150.0]))
# We have no grid
ws.VectorSet(ws.lat_grid, np.array([0.0]))
# We have no grid
ws.VectorSet(ws.lon_grid, np.array([0.0]))
# We have no grid
ws.IndexSet(ws.atmosphere_dim, 3)
# We have no atmosphere
ws.MatrixSet(ws.sensor_pos, np.array([[0.0, 0.0, 0.0]]))
# We have no sensor
ws.sensorOff()
# We have no sensor
ws.IndexSet(ws.propmat_clearsky_agenda_checked, 1)
# We have no propmat agenda
## Set up partial derivatives
ws.jacobianInit()
ws.jacobianAddTemperature(g1=np.array([150.0]), g2=np.array([0.0]), g3=np.array([0.0]))
ws.jacobianAddAbsSpecies(
    g1=np.array([150.0]),
    g2=np.array([0.0]),
    g3=np.array([0.0]),
    species="O2-66",
    for_species_tag=0,
)
ws.jacobianAddWind(
    g1=np.array([150.0]), g2=np.array([0.0]), g3=np.array([0.0]), dfrequency=0.1
)
ws.jacobianAddMagField(
    g1=np.array([150.0]), g2=np.array([0.0]), g3=np.array([0.0]), component="u"
)
ws.jacobianAddMagField(
    g1=np.array([150.0]), g2=np.array([0.0]), g3=np.array([0.0]), component="v"
)
ws.jacobianAddMagField(
    g1=np.array([150.0]), g2=np.array([0.0]), g3=np.array([0.0]), component="w"
)
ws.jacobianAddMagField(
    g1=np.array([150.0]), g2=np.array([0.0]), g3=np.array([0.0]), component="strength"
)
ws.jacobianClose()
ws.propmat_clearskyInit()
ws.lbl_checkedCalc()
ws.propmat_clearskyAddZeeman()
# WriteXML("ascii", propmat_clearsky, "testdata/zeeman/propmat.xml")
# WriteXML("ascii", dpropmat_clearsky_dx, "testdata/zeeman/dpropmat.xml")
ws.ReadXML(ws.test, "testdata/zeeman/propmat.xml")
ws.CompareRelative(ws.test, ws.propmat_clearsky, 1e-06)
ws.ReadXML(ws.test, "testdata/zeeman/dpropmat.xml")
ws.CompareRelative(ws.test, ws.dpropmat_clearsky_dx, 1e-06)
ws.jacobianOff()
ws.NumericSet(ws.rtp_temperature, 215.01)
ws.propmat_clearskyInit()
ws.propmat_clearskyAddZeeman()
# WriteXML("ascii", propmat_clearsky, "testdata/zeeman/propmat_dT.xml")
ws.ReadXML(ws.test, "testdata/zeeman/propmat_dT.xml")
ws.CompareRelative(ws.test, ws.propmat_clearsky, 1e-06)
ws.NumericSet(ws.rtp_temperature, 215.0)
ws.VectorSet(ws.rtp_vmr, np.array([0.21001]))
ws.propmat_clearskyInit()
ws.propmat_clearskyAddZeeman()
# WriteXML("ascii", propmat_clearsky, "testdata/zeeman/propmat_dvmr.xml")
ws.ReadXML(ws.test, "testdata/zeeman/propmat_dvmr.xml")
ws.CompareRelative(ws.test, ws.propmat_clearsky, 1e-06)
ws.VectorSet(ws.rtp_vmr, np.array([0.21]))
ws.VectorNLinSpace(ws.f_grid, 501, 99990100000.0, 100010100000.0)
ws.propmat_clearskyInit()
ws.propmat_clearskyAddZeeman()
# WriteXML("ascii", propmat_clearsky, "testdata/zeeman/propmat_df.xml")
ws.ReadXML(ws.test, "testdata/zeeman/propmat_df.xml")
ws.CompareRelative(ws.test, ws.propmat_clearsky, 1e-06)
ws.VectorNLinSpace(ws.f_grid, 501, 99990000000.0, 100010000000.0)
ws.VectorSet(ws.rtp_mag, np.array([2.501e-05, 6.000e-05, 1.000e-05]))
ws.propmat_clearskyInit()
ws.propmat_clearskyAddZeeman()
# WriteXML("ascii", propmat_clearsky, "testdata/zeeman/propmat_du.xml")
ws.ReadXML(ws.test, "testdata/zeeman/propmat_du.xml")
ws.CompareRelative(ws.test, ws.propmat_clearsky, 1e-06)
ws.VectorSet(ws.rtp_mag, np.array([2.500e-05, 6.001e-05, 1.000e-05]))
ws.propmat_clearskyInit()
ws.propmat_clearskyAddZeeman()
# WriteXML("ascii", propmat_clearsky, "testdata/zeeman/propmat_dv.xml")
ws.ReadXML(ws.test, "testdata/zeeman/propmat_dv.xml")
ws.CompareRelative(ws.test, ws.propmat_clearsky, 1e-06)
ws.VectorSet(ws.rtp_mag, np.array([2.500e-05, 6.000e-05, 1.001e-05]))
ws.propmat_clearskyInit()
ws.propmat_clearskyAddZeeman()
# WriteXML("ascii", propmat_clearsky, "testdata/zeeman/propmat_dw.xml")
ws.ReadXML(ws.test, "testdata/zeeman/propmat_dw.xml")
ws.CompareRelative(ws.test, ws.propmat_clearsky, 1e-06)
ws.VectorSet(ws.rtp_mag, np.array([2.5e-05, 6.0e-05, 1.0e-05]))
ws.VectorScale(ws.rtp_mag, ws.rtp_mag, 1.001)
ws.propmat_clearskyInit()
ws.propmat_clearskyAddZeeman()
# WriteXML("ascii", propmat_clearsky, "testdata/zeeman/propmat_dH.xml")
ws.ReadXML(ws.test, "testdata/zeeman/propmat_dH.xml")
ws.CompareRelative(ws.test, ws.propmat_clearsky, 1e-06)
ws.VectorSet(ws.rtp_mag, np.array([2.5e-05, 6.0e-05, 1.0e-05]))
