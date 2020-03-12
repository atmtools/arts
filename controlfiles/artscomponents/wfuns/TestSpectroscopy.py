import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)


@arts_agenda
def abs_xsec_agenda(ws):
    ws.abs_xsec_per_speciesInit()
    ws.abs_xsec_per_speciesAddLines()


ws.abs_xsec_agenda = abs_xsec_agenda


@arts_agenda
def propmat_clearsky_agenda(ws):
    ws.Ignore(ws.rtp_mag)
    ws.Ignore(ws.rtp_los)
    ws.propmat_clearskyInit()
    ws.propmat_clearskyAddOnTheFly()


ws.propmat_clearsky_agenda = propmat_clearsky_agenda


@arts_agenda
def iy_main_agenda(ws):
    ws.ppathCalc()
    ws.iyEmissionStandard()


ws.iy_main_agenda = iy_main_agenda


@arts_agenda
def geo_pos_agenda(ws):
    ws.Ignore(ws.ppath)
    ws.VectorSet(ws.geo_pos, [])


ws.geo_pos_agenda = geo_pos_agenda


@arts_agenda
def ppath_agenda(ws):
    ws.Ignore(ws.rte_pos2)
    ws.ppathStepByStep()


ws.ppath_agenda = ppath_agenda


@arts_agenda
def ppath_step_agenda(ws):
    ws.Ignore(ws.f_grid)
    ws.Ignore(ws.ppath_lraytrace)
    ws.ppath_stepGeometric()


ws.ppath_step_agenda = ppath_step_agenda


@arts_agenda
def water_p_eq_agenda(ws):
    ws.water_p_eq_fieldMK05()


ws.water_p_eq_agenda = water_p_eq_agenda


@arts_agenda
def iy_space_agenda(ws):
    ws.Ignore(ws.rtp_pos)
    ws.Ignore(ws.rtp_los)
    ws.MatrixCBR(ws.iy, ws.stokes_dim, ws.f_grid)


ws.iy_space_agenda = iy_space_agenda


@arts_agenda
def iy_surface_agenda(ws):
    ws.SurfaceDummy()
    ws.iySurfaceRtpropAgenda()


ws.iy_surface_agenda = iy_surface_agenda


@arts_agenda
def surface_rtprop_agenda(ws):
    ws.InterpAtmFieldToPosition(out=ws.surface_skin_t, field=ws.t_field)
    ws.surfaceBlackbody()


ws.surface_rtprop_agenda = surface_rtprop_agenda

ws.refellipsoidEarth()
## Constants (changes here might require changes in perturbations to match)
ws.isotopologue_ratiosInitFromBuiltin()
ws.ReadXML(ws.partition_functions, "spectroscopy/PartitionSums/TIPS/tips.xml")
ws.abs_speciesSet(species=["O2-66"])
ws.VectorNLinSpace(ws.f_grid, 100, 85000000000.0, 115000000000.0)
# FIXME: 101 freqs fail because HTP fails
ws.VectorNLogSpace(ws.p_grid, 51, 103000.0, 1.0)
ws.VectorSet(ws.lat_grid, np.array([-1.0, 1.0]))
# We have no grid
ws.VectorSet(ws.lon_grid, np.array([-1.0, 1.0]))
# We have no grid
ws.z_surfaceConstantAltitude()
ws.NumericSet(ws.rte_alonglos_v, 0.0)
ws.NumericSet(ws.ppath_lmax, 1000.0)
ws.NumericSet(ws.ppath_lraytrace, 1000.0)
ws.IndexSet(ws.stokes_dim, 1)
ws.IndexSet(ws.abs_f_interp_order, 1)
ws.Touch(ws.wind_u_field)
ws.Touch(ws.wind_v_field)
ws.Touch(ws.wind_w_field)
ws.Touch(ws.mag_u_field)
ws.Touch(ws.mag_v_field)
ws.Touch(ws.mag_w_field)
ws.Touch(ws.iy_aux_vars)
ws.Touch(ws.surface_props_data)
ws.Touch(ws.surface_props_names)
ws.Touch(ws.transmitter_pos)
ws.StringSet(ws.iy_unit, "PlanckBT")
ws.AtmRawRead(basename="testdata/tropical")
ws.AtmosphereSet3D()
ws.AtmFieldsCalcExpand1D()
## Calculate w/o NLTE
ws.nlteOff()
## Comparative parameter
ws.VectorCreate("testy")
ws.MatrixCreate("testjac")
## Absorption lines (Doppler, no Zeeman, no NLTE, no mirroring, no normalization)
ws.ReadXML(ws.abs_lines, "../lineshapes/testdata/lm-htp-line.xml")
ws.ArrayOfAbsorptionLinesCreate("aolr")
ws.Copy(ws.aolr, ws.abs_lines)
ws.abs_lines_per_speciesCreateFromLines()
## Line matching information
ws.ArrayOfQuantumIdentifierCreate("qi_lines")
ws.ReadXML(ws.qi_lines, "../lineshapes/testdata/qi-line.xml")
ws.QuantumIdentifierCreate("QI")
ws.Extract(ws.QI, ws.qi_lines, 0)
ws.MatrixSet(ws.sensor_pos, np.array([[300000.0, 0.0, 0.0]]))
# 300 km altitude
ws.MatrixSet(ws.sensor_los, np.array([[180.0, 0.0]]))
# Nadir looking
ws.sensorOff()
# We have no sensor
## Set up partial derivatives
ws.jacobianInit()
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
ws.jacobianAddShapeCatalogParameter(
    line_identity=ws.QI, variable="G2", coefficient="X0", species="AIR"
)
ws.jacobianAddShapeCatalogParameter(
    line_identity=ws.QI, variable="G2", coefficient="X1", species="AIR"
)
ws.jacobianAddShapeCatalogParameter(
    line_identity=ws.QI, variable="D2", coefficient="X0", species="AIR"
)
ws.jacobianAddShapeCatalogParameter(
    line_identity=ws.QI, variable="D2", coefficient="X1", species="AIR"
)
ws.jacobianAddShapeCatalogParameter(
    line_identity=ws.QI, variable="FVC", coefficient="X0", species="AIR"
)
ws.jacobianAddShapeCatalogParameter(
    line_identity=ws.QI, variable="FVC", coefficient="X1", species="AIR"
)
ws.jacobianAddShapeCatalogParameter(
    line_identity=ws.QI, variable="ETA", coefficient="X0", species="AIR"
)
ws.jacobianAddShapeCatalogParameter(
    line_identity=ws.QI, variable="ETA", coefficient="X1", species="AIR"
)
ws.jacobianAddShapeCatalogParameter(
    line_identity=ws.QI, variable="Y", coefficient="X0", species="AIR"
)
ws.jacobianAddShapeCatalogParameter(
    line_identity=ws.QI, variable="Y", coefficient="X1", species="AIR"
)
ws.jacobianAddShapeCatalogParameter(
    line_identity=ws.QI, variable="Y", coefficient="X2", species="AIR"
)
ws.jacobianAddShapeCatalogParameter(
    line_identity=ws.QI, variable="G", coefficient="X0", species="AIR"
)
ws.jacobianAddShapeCatalogParameter(
    line_identity=ws.QI, variable="G", coefficient="X1", species="AIR"
)
ws.jacobianAddShapeCatalogParameter(
    line_identity=ws.QI, variable="G", coefficient="X2", species="AIR"
)
ws.jacobianAddShapeCatalogParameter(
    line_identity=ws.QI, variable="DV", coefficient="X0", species="AIR"
)
ws.jacobianAddShapeCatalogParameter(
    line_identity=ws.QI, variable="DV", coefficient="X1", species="AIR"
)
ws.jacobianAddShapeCatalogParameter(
    line_identity=ws.QI, variable="DV", coefficient="X2", species="AIR"
)
ws.jacobianClose()
ws.cloudboxOff()
# Perform calculations for analytical propagation matrix and derivatives
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
ws.sensor_checkedCalc()
ws.lbl_checkedCalc()
ws.yCalc()
# WriteXML("ascii", y, "testdata/comparedata/y.xml")
# WriteXML("ascii", jacobian, "testdata/comparedata/dy.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
ws.ReadXML(ws.testjac, "testdata/comparedata/dy.xml")
ws.CompareRelative(ws.testjac, ws.jacobian, 0.2)
# Turn off the jacobian to make for faster calculations for perturbations below
ws.jacobianOff()
ws.cloudboxOff()
# Perform calculations for perturbed line strength derivative
ws.abs_linesChangeBaseParameterForMatchingLines(
    QI=ws.QI, parameter_name="Line Strength", change=1e-17, relative=0, loose_matching=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-ds0.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-ds0.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
# Perform calculations for perturbed line center derivative
ws.abs_linesChangeBaseParameterForMatchingLines(
    QI=ws.QI, parameter_name="Line Center", change=10.0, relative=0, loose_matching=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-df0.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-df0.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
# Perform calculations for perturbed AIR-G0-X0 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="G0", coefficient="X0", species="AIR", change=20.0, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-dAIR-G0-X0.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-dAIR-G0-X0.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
# Perform calculations for perturbed AIR-G0-X1 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="G0", coefficient="X1", species="AIR", change=0.0008, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-dAIR-G0-X1.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-dAIR-G0-X1.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
# Perform calculations for perturbed AIR-D0-X0 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="D0", coefficient="X0", species="AIR", change=1.0, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-dAIR-D0-X0.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-dAIR-D0-X0.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
# Perform calculations for perturbed AIR-D0-X1 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="D0", coefficient="X1", species="AIR", change=0.0008, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-dAIR-D0-X1.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-dAIR-D0-X1.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
# Perform calculations for perturbed AIR-G2-X0 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="G2", coefficient="X0", species="AIR", change=4.0, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-dAIR-G2-X0.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-dAIR-G2-X0.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
# Perform calculations for perturbed AIR-G2-X1 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="G2", coefficient="X1", species="AIR", change=0.01, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-dAIR-G2-X1.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-dAIR-G2-X1.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
# Perform calculations for perturbed AIR-D2-X0 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="D2", coefficient="X0", species="AIR", change=0.5, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-dAIR-D2-X0.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-dAIR-D2-X0.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
# Perform calculations for perturbed AIR-D2-X1 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="D2", coefficient="X1", species="AIR", change=0.01, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-dAIR-D2-X1.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-dAIR-D2-X1.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
# Perform calculations for perturbed AIR-FVC-X0 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="FVC", coefficient="X0", species="AIR", change=20.0, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-dAIR-FVC-X0.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-dAIR-FVC-X0.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
# Perform calculations for perturbed AIR-FVC-X1 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="FVC", coefficient="X1", species="AIR", change=2.0, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-dAIR-FVC-X1.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-dAIR-FVC-X1.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
# Perform calculations for perturbed AIR-ETA-X0 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI,
    parameter="ETA",
    coefficient="X0",
    species="AIR",
    change=0.0001,
    relative=0,
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-dAIR-ETA-X0.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-dAIR-ETA-X0.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
# Perform calculations for perturbed AIR-ETA-X1 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI,
    parameter="ETA",
    coefficient="X1",
    species="AIR",
    change=0.0001,
    relative=0,
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-dAIR-ETA-X1.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-dAIR-ETA-X1.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
# Perform calculations for perturbed AIR-ETA-X0 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="Y", coefficient="X0", species="AIR", change=1e-11, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-dAIR-Y-X0.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-dAIR-Y-X0.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
# Perform calculations for perturbed AIR-ETA-X0 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="Y", coefficient="X1", species="AIR", change=1e-13, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-dAIR-Y-X1.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-dAIR-Y-X1.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
# Perform calculations for perturbed AIR-ETA-X0 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="Y", coefficient="X2", species="AIR", change=0.01, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-dAIR-Y-X2.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-dAIR-Y-X2.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
# Perform calculations for perturbed AIR-ETA-X0 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="G", coefficient="X0", species="AIR", change=1e-15, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-dAIR-G-X0.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-dAIR-G-X0.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
# Perform calculations for perturbed AIR-ETA-X0 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="G", coefficient="X1", species="AIR", change=1e-17, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-dAIR-G-X1.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-dAIR-G-X1.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
# Perform calculations for perturbed AIR-ETA-X0 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="G", coefficient="X2", species="AIR", change=0.01, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-dAIR-G-X2.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-dAIR-G-X2.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
ws.CompareRelative(ws.testy, ws.y, 1e-06)
# Perform calculations for perturbed AIR-ETA-X0 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="DV", coefficient="X0", species="AIR", change=1e-07, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-dAIR-DV-X0.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-dAIR-DV-X0.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
# Perform calculations for perturbed AIR-ETA-X0 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="DV", coefficient="X1", species="AIR", change=1e-09, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-dAIR-DV-X1.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-dAIR-DV-X1.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
# Perform calculations for perturbed AIR-ETA-X0 derivative
ws.abs_linesChangeLineShapeModelParameterForMatchingLines(
    QI=ws.QI, parameter="DV", coefficient="X2", species="AIR", change=0.01, relative=0
)
ws.abs_lines_per_speciesCreateFromLines()
ws.abs_xsec_agenda_checkedCalc()
ws.yCalc()
ws.Copy(ws.abs_lines, ws.aolr)
# WriteXML("ascii", y, "testdata/comparedata/y-dAIR-DV-X2.xml")
ws.ReadXML(ws.testy, "testdata/comparedata/y-dAIR-DV-X2.xml")
ws.CompareRelative(ws.testy, ws.y, 1e-06)
