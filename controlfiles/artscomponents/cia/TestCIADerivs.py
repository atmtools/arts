import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)


@arts_agenda
def abs_xsec_agenda(ws):
    ws.abs_xsec_per_speciesInit()
    ws.abs_xsec_per_speciesAddCIA()
    ws.abs_xsec_per_speciesAddLines()


ws.abs_xsec_agenda = abs_xsec_agenda

#   Can run with abs_xsec_per_speciesAddLines by uncommenting:
#   abs_lineshapeDefine(abs_lineshape, "Faddeeva_Algorithm_916", "no_norm", -1)
#   Will fail at first derivatives with somewhat small error
# Constants (changes here might require changes in perturbations to match)
ws.isotopologue_ratiosInitFromBuiltin()
ws.partition_functionsInitFromBuiltin()
ws.abs_speciesSet(species=["H2-CIA-H2-0,H2-CIA-He-0", "He"])
# BUG:  Must have lines (though good for this test)
ws.VectorNLinSpace(ws.f_grid, 11, 3000000000.0, 30000000000000.0)
ws.VectorSet(ws.rtp_vmr, array([0.1, 0.2]))
ws.Touch(ws.rtp_nlte)
ws.NumericSet(ws.rtp_temperature, 250.0)
ws.NumericSet(ws.rtp_pressure, 25000.0)
ws.IndexSet(ws.stokes_dim, 1)
ws.Touch(ws.rtp_nlte)
ws.IndexSet(ws.abs_f_interp_order, 3)
# Calculate w/o NLTE
ws.nlteOff()
# Comparative parameter
ws.ArrayOfPropagationMatrixCreate("testdata")
# CIA catalog from arts-xml-data
ws.ReadXML(ws.abs_cia_data, "spectroscopy/cia/borysow/Borysow_CIA_JUICE_SWI.xml")
# Absorption lines (VP)
ws.ReadARTSCAT(ws.abs_lines, "tests/test-he-line.xml")
ws.abs_lines_per_speciesCreateFromLines()
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
ws.jacobianAddAbsSpecies(
    g1=ws.p_grid, g2=array([0.0]), g3=array([0.0]), species="H2", for_species_tag=0
)
ws.jacobianAddAbsSpecies(
    g1=ws.p_grid, g2=array([0.0]), g3=array([0.0]), species="He", for_species_tag=0
)
ws.jacobianAddWind(g1=ws.p_grid, g2=array([0.0]), g3=array([0.0]), dfrequency=1000.0)
ws.jacobianClose()
# Perform calculations for analytical propagation matrix and derivatives
ws.abs_xsec_agenda_checkedCalc()
ws.lbl_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
#   WriteXML("ascii", propmat_clearsky, "tests/propmat.xml")
#   WriteXML("ascii", dpropmat_clearsky_dx, "tests/dpropmat.xml")
ws.ReadXML(ws.testdata, "tests/propmat.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
ws.ReadXML(ws.testdata, "tests/dpropmat.xml")
ws.CompareRelative(ws.testdata, ws.dpropmat_clearsky_dx, 0.0001)
# Turn off the jacobian to make for faster calculations for perturbations below
ws.jacobianOff()
# Perform calculations for perturbed temperature derivative
ws.NumericSet(ws.rtp_temperature, 250.0001)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.NumericSet(ws.rtp_temperature, 250.0)
#   WriteXML("ascii", propmat_clearsky, "tests/propmat-dT.xml")
ws.ReadXML(ws.testdata, "tests/propmat-dT.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed VMR derivative
ws.VectorSet(ws.rtp_vmr, array([0.10001, 0.2]))
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.VectorSet(ws.rtp_vmr, array([0.1, 0.2]))
#   WriteXML("ascii", propmat_clearsky, "tests/propmat-dvmr-H2.xml")
ws.ReadXML(ws.testdata, "tests/propmat-dvmr-H2.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed VMR derivative
ws.VectorSet(ws.rtp_vmr, array([0.1, 0.20001]))
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.VectorSet(ws.rtp_vmr, array([0.1, 0.2]))
#   WriteXML("ascii", propmat_clearsky, "tests/propmat-dvmr-He.xml")
ws.ReadXML(ws.testdata, "tests/propmat-dvmr-He.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
# Perform calculations for perturbed frequency derivative
ws.VectorNLinSpace(ws.f_grid, 11, 3000001000.0, 30000000001000.0)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearskyInit()
ws.propmat_clearskyAddOnTheFly()
ws.VectorNLinSpace(ws.f_grid, 11, 3000000000.0, 30000000000000.0)
#   WriteXML("ascii", propmat_clearsky, "tests/propmat-df.xml")
ws.ReadXML(ws.testdata, "tests/propmat-df.xml")
ws.CompareRelative(ws.testdata, ws.propmat_clearsky, 1e-06)
