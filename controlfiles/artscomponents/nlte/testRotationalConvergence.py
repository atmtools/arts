import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)


@arts_agenda
def iy_space_agenda(ws):
    ws.Ignore(ws.rtp_pos)
    ws.Ignore(ws.rtp_los)
    ws.MatrixCBR(ws.iy, ws.stokes_dim, ws.f_grid)


ws.iy_space_agenda = iy_space_agenda


@arts_agenda
def iy_main_agenda(ws):
    ws.ppathCalc()
    ws.iyEmissionStandard()


ws.iy_main_agenda = iy_main_agenda


@arts_agenda
def water_p_eq_agenda(ws):
    ws.water_p_eq_fieldMK05()


ws.water_p_eq_agenda = water_p_eq_agenda


@arts_agenda
def iy_surface_agenda(ws):
    ws.SurfaceDummy()
    ws.iySurfaceRtpropAgenda()


ws.iy_surface_agenda = iy_surface_agenda


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


@arts_agenda
def surface_rtprop_agenda(ws):
    ws.InterpAtmFieldToPosition(out=ws.surface_skin_t, field=ws.t_field)
    ws.surfaceBlackbody()


ws.surface_rtprop_agenda = surface_rtprop_agenda

# Names that are uninitialized but needs to be initialized empty for the code to work
ws.Touch(ws.surface_props_names)
ws.Touch(ws.surface_props_data)
ws.Touch(ws.wind_u_field)
ws.Touch(ws.wind_v_field)
ws.Touch(ws.wind_w_field)
ws.Touch(ws.mag_u_field)
ws.Touch(ws.mag_v_field)
ws.Touch(ws.mag_w_field)
# Set constants
ws.IndexSet(ws.stokes_dim, 1)
ws.NumericSet(ws.rte_alonglos_v, 0.0)
ws.NumericSet(ws.ppath_lmax, -1.0)
ws.MatrixSetConstant(ws.z_surface, 1, 1, 0.0)
ws.abs_speciesSet(species=["H2O-161"])
# Data reading
ws.AtmRawRead(basename="testdata/")
ws.ReadARTSCAT(
    abs_lines=ws.abs_lines,
    filename="testdata/lines.xml",
    globalquantumnumbers="Ka Kc J",
)
ws.collision_coefficientsFromSplitFiles(basename="spectroscopy/nlte/H2O/")
ws.ReadXML(ws.nlte_level_identifiers, "testdata/qi.xml")
# Interpretation of data
ws.isotopologue_ratiosInitFromBuiltin()
ws.partition_functionsInitFromBuiltin()
ws.p_gridFromZRaw()
ws.AtmosphereSet1D()
ws.refellipsoidGanymede(model="Sphere")
ws.AtmFieldsCalc()
ws.abs_lines_per_speciesCreateFromLines()
# No fancy computations by turning all of these off
ws.jacobianOff()
ws.cloudboxOff()
# Initialize the ratio-field
ws.nlte_fieldSetLteInternalPartitionFunction()
ws.nlte_fieldRescalePopulationLevels(s=0.75)
# Check that we are performing sane computations
ws.propmat_clearsky_agenda_checkedCalc()
ws.abs_xsec_agenda_checkedCalc()
ws.lbl_checkedCalc()
# Compute the new nlte_field
ws.nlte_fieldForSingleSpeciesNonOverlappingLines(
    df=0.0001, nz=10, nf=401, dampened=0, iteration_limit=100, convergence_limit=0.0001
)
# WriteXML("ascii", nlte_field, "testdata/nlte_testdata.xml")
# Compare with stored data to make sure it is OK
ws.EnergyLevelMapCreate("nlte_testdata")
ws.ReadXML(ws.nlte_testdata, "testdata/nlte_testdata.xml")
ws.CompareRelative(ws.nlte_testdata, ws.nlte_field, 1e-06)
