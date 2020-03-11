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
# Constants (changes here might require changes in perturbations to match)
ws.isotopologue_ratiosInitFromBuiltin()
ws.partition_functionsInitFromBuiltin()
ws.abs_speciesSet(species=["O2-66"])
ws.VectorNLinSpace(ws.f_grid, 101, 99990000000.0, 100010000000.0)
ws.VectorSet(ws.rtp_vmr, array([0.21]))
ws.NumericSet(ws.rtp_temperature, 250.0)
ws.NumericSet(ws.rtp_pressure, 1e-05)
ws.IndexSet(ws.stokes_dim, 1)
# Comparative parameter
ws.ArrayOfPropagationMatrixCreate("testdata")
# Absorption line (Doppler, no Zeeman, no NLTE, no mirroring, no normalization)
ws.ReadXML(ws.abs_lines, "../lineshapes/testdata/lm-vp-line.xml")
ws.ArrayOfAbsorptionLinesCreate("aolr")
ws.Copy(ws.aolr, ws.abs_lines)
ws.abs_lines_per_speciesCreateFromLines()
# Line matching information
ws.ArrayOfQuantumIdentifierCreate("qi_lines")
ws.ReadXML(ws.qi_lines, "../lineshapes/testdata/qi-line.xml")
ws.QuantumIdentifierCreate("QI")
ws.Extract(ws.QI, ws.qi_lines, 0)
# NLTE matching information
ws.ReadXML(ws.nlte_level_identifiers, "testdata/o2-qi.xml")
ws.Touch(ws.nlte_vibrational_energies)
ws.rtp_nlteFromRaw(data=array([0.7, 0.4]))
ws.nlteSetByQuantumIdentifiers(nlte_field=ws.rtp_nlte)
# Silly parameters that have to be set by agendas and ARTS in general but are completely useless for these calculations
ws.VectorSet(ws.p_grid, array([150.0]))
# We have no grid
ws.VectorSet(ws.lat_grid, array([0.0]))
# We have no grid
ws.VectorSet(ws.lon_grid, array([0.0]))
# We have no grid
ws.IndexSet(ws.atmosphere_dim, 1)
# We have no atmosphere
ws.MatrixSet(ws.sensor_pos, array([[0.0, 0.0, 0.0]]))
# We have no sensor
ws.sensorOff()
# We have no sensor
ws.IndexSet(ws.propmat_clearsky_agenda_checked, 1)
# We have no propmat agenda
# Set up partial derivatives
ws.jacobianInit()
ws.jacobianAddTemperature(g1=ws.p_grid, g2=array([0.0]), g3=array([0.0]))
ws.jacobianAddWind(g1=ws.p_grid, g2=array([0.0]), g3=array([0.0]), dfrequency=0.1)
ws.jacobianAddAbsSpecies(
    g1=ws.p_grid, g2=array([0.0]), g3=array([0.0]), species="O2-66", for_species_tag=0
)
ws.jacobianAddBasicCatalogParameters(
    catalog_identities=ws.qi_lines, catalog_parameters=["Line Strength", "Line Center"]
)
ws.jacobianAddShapeCatalogParameters(
    line_identities=ws.qi_lines,
    variables=["G0", "D0", "Y", "G", "DV"],
    coefficients=["X0", "X1", "X2"],
    species=["AIR"],
)
ws.jacobianAddNLTEs(
    g1=ws.p_grid,
    g2=array([0.0]),
    g3=array([0.0]),
    energy_level_identities=ws.nlte_level_identifiers,
    dx=1e-06,
)
ws.jacobianClose()
# Perform calculations for analytical propagation matrix and derivatives
ws.abs_xsec_agenda_checkedCalc()
ws.lbl_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.lbl_checkedCalc()
# WriteXML("ascii", propmat_clearsky, "testdata/test-nlte/propmat.xml")
# WriteXML("ascii", dpropmat_clearsky_dx, "testdata/test-nlte/dpropmat.xml")
ws.ReadXML(ws.testdata, "testdata/test-nlte/propmat.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
ws.ReadXML(ws.testdata, "testdata/test-nlte/dpropmat.xml")
ws.CompareRelative(ws.testdata, ws.dpropmat_clearsky_dx, 0.0001)
# Turn off the jacobian to make for faster calculations for perturbations below
ws.jacobianOff()
# Perform calculations for perturbed temperature derivative
ws.NumericSet(ws.rtp_temperature, 250.0001)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.NumericSet(ws.rtp_temperature, 250.0)
# WriteXML("ascii", propmat_clearsky, "testdata/test-nlte/propmat-dT.xml")
ws.ReadXML(ws.testdata, "testdata/test-nlte/propmat-dT.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed VMR derivative
ws.VectorSet(ws.rtp_vmr, array([0.2101]))
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.VectorSet(ws.rtp_vmr, array([0.21]))
# WriteXML("ascii", propmat_clearsky, "testdata/test-nlte/propmat-dvmr.xml")
ws.ReadXML(ws.testdata, "testdata/test-nlte/propmat-dvmr.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed frequency derivative
ws.VectorNLinSpace(ws.f_grid, 101, 99990000100.0, 100010000100.0)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.VectorNLinSpace(ws.f_grid, 101, 99990000000.0, 100010000000.0)
# WriteXML("ascii", propmat_clearsky, "testdata/test-nlte/propmat-df.xml")
ws.ReadXML(ws.testdata, "testdata/test-nlte/propmat-df.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed line strength derivative
ws.abs_linesChangeBaseParameterForMatchingLines(
    QI=ws.QI, parameter_name="Line Strength", change=1e-30, relative=0, loose_matching=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.nlteSetByQuantumIdentifiers(nlte_field=ws.rtp_nlte)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-nlte/propmat-ds0.xml")
ws.ReadXML(ws.testdata, "testdata/test-nlte/propmat-ds0.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed line center derivative
ws.abs_linesChangeBaseParameterForMatchingLines(
    QI=ws.QI, parameter_name="Line Center", change=10.0, relative=0, loose_matching=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.nlteSetByQuantumIdentifiers(nlte_field=ws.rtp_nlte)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-nlte/propmat-df0.xml")
ws.ReadXML(ws.testdata, "testdata/test-nlte/propmat-df0.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed AIR-G0-X0 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="G0", coefficient="X0", species="AIR", change=10.0, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.nlteSetByQuantumIdentifiers(nlte_field=ws.rtp_nlte)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-nlte/propmat-dlf-AIR-G0-X0.xml")
ws.ReadXML(ws.testdata, "testdata/test-nlte/propmat-dlf-AIR-G0-X0.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed AIR-G0-X1 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="G0", coefficient="X1", species="AIR", change=1e-05, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.nlteSetByQuantumIdentifiers(nlte_field=ws.rtp_nlte)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-nlte/propmat-dlf-AIR-G0-X1.xml")
ws.ReadXML(ws.testdata, "testdata/test-nlte/propmat-dlf-AIR-G0-X1.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed AIR-D0-X0 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="D0", coefficient="X0", species="AIR", change=10.0, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.nlteSetByQuantumIdentifiers(nlte_field=ws.rtp_nlte)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-nlte/propmat-dlf-AIR-D0-X0.xml")
ws.ReadXML(ws.testdata, "testdata/test-nlte/propmat-dlf-AIR-D0-X0.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed AIR-D0-X1 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="D0", coefficient="X1", species="AIR", change=1e-05, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.nlteSetByQuantumIdentifiers(nlte_field=ws.rtp_nlte)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-nlte/propmat-dlf-AIR-D0-X1.xml")
ws.ReadXML(ws.testdata, "testdata/test-nlte/propmat-dlf-AIR-D0-X1.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed AIR-Y-X0 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="Y", coefficient="X0", species="AIR", change=1e-10, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.nlteSetByQuantumIdentifiers(nlte_field=ws.rtp_nlte)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-nlte/propmat-dlf-AIR-Y-X0.xml")
ws.ReadXML(ws.testdata, "testdata/test-nlte/propmat-dlf-AIR-Y-X0.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed AIR-Y-X1 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="Y", coefficient="X1", species="AIR", change=1e-12, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.nlteSetByQuantumIdentifiers(nlte_field=ws.rtp_nlte)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-nlte/propmat-dlf-AIR-Y-X1.xml")
ws.ReadXML(ws.testdata, "testdata/test-nlte/propmat-dlf-AIR-Y-X1.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed AIR-Y-X2 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="Y", coefficient="X2", species="AIR", change=1e-05, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.nlteSetByQuantumIdentifiers(nlte_field=ws.rtp_nlte)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-nlte/propmat-dlf-AIR-Y-X2.xml")
ws.ReadXML(ws.testdata, "testdata/test-nlte/propmat-dlf-AIR-Y-X2.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed AIR-G-X0 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="G", coefficient="X0", species="AIR", change=1e-14, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.nlteSetByQuantumIdentifiers(nlte_field=ws.rtp_nlte)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-nlte/propmat-dlf-AIR-G-X0.xml")
ws.ReadXML(ws.testdata, "testdata/test-nlte/propmat-dlf-AIR-G-X0.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed AIR-G-X1 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="G", coefficient="X1", species="AIR", change=1e-16, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.nlteSetByQuantumIdentifiers(nlte_field=ws.rtp_nlte)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-nlte/propmat-dlf-AIR-G-X1.xml")
ws.ReadXML(ws.testdata, "testdata/test-nlte/propmat-dlf-AIR-G-X1.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed AIR-G-X2 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="G", coefficient="X2", species="AIR", change=1e-05, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.nlteSetByQuantumIdentifiers(nlte_field=ws.rtp_nlte)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-nlte/propmat-dlf-AIR-G-X2.xml")
ws.ReadXML(ws.testdata, "testdata/test-nlte/propmat-dlf-AIR-G-X2.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed AIR-DV-X0 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="DV", coefficient="X0", species="AIR", change=0.1, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.nlteSetByQuantumIdentifiers(nlte_field=ws.rtp_nlte)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-nlte/propmat-dlf-AIR-DV-X0.xml")
ws.ReadXML(ws.testdata, "testdata/test-nlte/propmat-dlf-AIR-DV-X0.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed AIR-DV-X1 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="DV", coefficient="X1", species="AIR", change=0.01, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.nlteSetByQuantumIdentifiers(nlte_field=ws.rtp_nlte)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-nlte/propmat-dlf-AIR-DV-X1.xml")
ws.ReadXML(ws.testdata, "testdata/test-nlte/propmat-dlf-AIR-DV-X1.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed AIR-DV-X2 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="DV", coefficient="X2", species="AIR", change=1e-05, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.nlteSetByQuantumIdentifiers(nlte_field=ws.rtp_nlte)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-nlte/propmat-dlf-AIR-DV-X2.xml")
ws.ReadXML(ws.testdata, "testdata/test-nlte/propmat-dlf-AIR-DV-X2.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed lower levels derivative
ws.rtp_nlteFromRaw(data=array([0.7, 0.40001]))
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.rtp_nlteFromRaw(data=array([0.7, 0.4]))
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-nlte/propmat-dlow.xml")
ws.ReadXML(ws.testdata, "testdata/test-nlte/propmat-dlow.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed lower levels derivative
ws.rtp_nlteFromRaw(data=array([0.70001, 0.4]))
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.rtp_nlteFromRaw(data=array([0.7, 0.4]))
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", propmat_clearsky, "testdata/test-nlte/propmat-dupp.xml")
ws.ReadXML(ws.testdata, "testdata/test-nlte/propmat-dupp.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
