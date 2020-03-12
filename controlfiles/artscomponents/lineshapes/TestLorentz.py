import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)


@arts_agenda
def abs_xsec_agenda(ws):
    ws.abs_xsec_per_speciesInit()
    ws.abs_xsec_per_speciesAddLines()


ws.abs_xsec_agenda = abs_xsec_agenda

#   Can run with abs_xsec_per_speciesAddLines by uncommenting:
#   abs_lineshapeDefine(abs_lineshape, "Faddeeva_Algorithm_916", "no_norm", -1)
#   Will fail at first derivatives...
## Constants (changes here might require changes in perturbations to match)
ws.isotopologue_ratiosInitFromBuiltin()
ws.partition_functionsInitFromBuiltin()
ws.abs_speciesSet(species=["O2-66"])
ws.VectorNLinSpace(ws.f_grid, 101, 90000000000.0, 110000000000.0)
ws.Touch(ws.rtp_nlte)
ws.VectorSet(ws.rtp_vmr, np.array([0.21]))
ws.NumericSet(ws.rtp_temperature, 250.0)
ws.NumericSet(ws.rtp_pressure, 25000.0)
ws.IndexSet(ws.stokes_dim, 1)
## Calculate w/o NLTE
ws.nlteOff()
## Comparative parameter
ws.ArrayOfPropagationMatrixCreate("testdata")
## Absorption lines (Doppler, no Zeeman, no NLTE, no mirroring, no normalization)
ws.ReadXML(ws.abs_lines, "testdata/lp-line.xml")
ws.ArrayOfAbsorptionLinesCreate("aolr")
ws.Copy(ws.aolr, ws.abs_lines)
ws.abs_lines_per_speciesCreateFromLines()
## Line matching information
ws.ArrayOfQuantumIdentifierCreate("qi_lines")
ws.ReadXML(ws.qi_lines, "testdata/qi-line.xml")
ws.QuantumIdentifierCreate("QI")
ws.Extract(ws.QI, ws.qi_lines, 0)
## Silly parameters that have to be set by agendas and ARTS in general but are completely useless for these calculations
ws.VectorSet(ws.p_grid, np.array([150.0]))
# We have no grid
ws.VectorSet(ws.lat_grid, np.array([0.0]))
# We have no grid
ws.VectorSet(ws.lon_grid, np.array([0.0]))
# We have no grid
ws.IndexSet(ws.atmosphere_dim, 1)
# We have no atmosphere
ws.MatrixSet(ws.sensor_pos, np.array([[0.0, 0.0, 0.0]]))
# We have no sensor
ws.sensorOff()
# We have no sensor
ws.IndexSet(ws.propmat_clearsky_agenda_checked, 1)
# We have no propmat agenda
## Set up partial derivatives
ws.jacobianInit()
ws.jacobianAddTemperature(g1=ws.p_grid, g2=np.array([0.0]), g3=np.array([0.0]))
ws.jacobianAddWind(g1=ws.p_grid, g2=np.array([0.0]), g3=np.array([0.0]), dfrequency=0.1)
ws.jacobianAddAbsSpecies(
    g1=ws.p_grid,
    g2=np.array([0.0]),
    g3=np.array([0.0]),
    species="O2-66",
    for_species_tag=0,
)
ws.jacobianAddBasicCatalogParameters(
    catalog_identities=ws.qi_lines, catalog_parameters=["Line Strength", "Line Center"]
)
ws.jacobianAddShapeCatalogParameter(
    line_identity=ws.QI, variable="G0", coefficient="X0", species="AIR"
)
ws.jacobianAddShapeCatalogParameter(
    line_identity=ws.QI, variable="G0", coefficient="X1", species="AIR"
)
ws.jacobianAddShapeCatalogParameter(
    line_identity=ws.QI, variable="D0", coefficient="X0", species="AIR"
)
ws.jacobianAddShapeCatalogParameter(
    line_identity=ws.QI, variable="D0", coefficient="X1", species="AIR"
)
ws.jacobianClose()
# Perform calculations for analytical propagation matrix and derivatives
ws.abs_xsec_agenda_checkedCalc()
ws.lbl_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
# WriteXML("ascii", propmat_clearsky, "testdata/test-lp/propmat.xml")
# WriteXML("ascii", dpropmat_clearsky_dx, "testdata/test-lp/dpropmat.xml")
ws.ReadXML(ws.testdata, "testdata/test-lp/propmat.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
ws.ReadXML(ws.testdata, "testdata/test-lp/dpropmat.xml")
ws.CompareRelative(ws.testdata, ws.dpropmat_clearsky_dx, 0.0001)
# Turn off the jacobian to make for faster calculations for perturbations below
ws.jacobianOff()
# Perform calculations for perturbed temperature derivative
ws.NumericSet(ws.rtp_temperature, 250.0001)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.NumericSet(ws.rtp_temperature, 250.0)
# WriteXML("ascii", propmat_clearsky, "testdata/test-lp/propmat-dT.xml")
ws.ReadXML(ws.testdata, "testdata/test-lp/propmat-dT.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed frequency derivative
ws.VectorNLinSpace(ws.f_grid, 101, 90000000100.0, 110000000100.0)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.VectorNLinSpace(ws.f_grid, 101, 90000000000.0, 110000000000.0)
# WriteXML("ascii", propmat_clearsky, "testdata/test-lp/propmat-df.xml")
ws.ReadXML(ws.testdata, "testdata/test-lp/propmat-df.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed VMR derivative
ws.VectorSet(ws.rtp_vmr, np.array([0.2101]))
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.VectorSet(ws.rtp_vmr, np.array([0.21]))
# WriteXML("ascii", propmat_clearsky, "testdata/test-lp/propmat-dvmr.xml")
ws.ReadXML(ws.testdata, "testdata/test-lp/propmat-dvmr.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed line strength derivative
ws.abs_linesChangeBaseParameterForMatchingLines(
    QI=ws.QI, parameter_name="Line Strength", change=1e-17, relative=0, loose_matching=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-lp/propmat-ds0.xml")
ws.ReadXML(ws.testdata, "testdata/test-lp/propmat-ds0.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed line center derivative
ws.abs_linesChangeBaseParameterForMatchingLines(
    QI=ws.QI, parameter_name="Line Center", change=10.0, relative=0, loose_matching=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-lp/propmat-df0.xml")
ws.ReadXML(ws.testdata, "testdata/test-lp/propmat-df0.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed AIR-G0-X0 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="G0", coefficient="X0", species="AIR", change=20.0, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-lp/propmat-G0-X0.xml")
ws.ReadXML(ws.testdata, "testdata/test-lp/propmat-G0-X0.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed AIR-G0-X1 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="G0", coefficient="X1", species="AIR", change=0.0008, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-lp/propmat-G0-X1.xml")
ws.ReadXML(ws.testdata, "testdata/test-lp/propmat-G0-X1.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed AIR-D0-X0 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="D0", coefficient="X0", species="AIR", change=1.0, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-lp/propmat-D0-X0.xml")
ws.ReadXML(ws.testdata, "testdata/test-lp/propmat-D0-X0.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed AIR-D0-X1 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="D0", coefficient="X1", species="AIR", change=0.0008, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-lp/propmat-D0-X1.xml")
ws.ReadXML(ws.testdata, "testdata/test-lp/propmat-D0-X1.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
