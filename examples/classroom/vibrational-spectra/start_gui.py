import pyarts
import numpy as np

ws = pyarts.workspace.Workspace()

# Species setup and catalog reading
ws.abs_speciesSet(
    species=["HCl",
             "ClO",
             "CO",
             "N2O",
             "O3",
             "H2O",
             "CO2",
             "O2",
             "CH4"])
ws.abs_lines_per_speciesReadSpeciesSplitCatalog(basename="lines/")

# Speed up by limiting the frequency range of lines
ws.abs_lines_per_speciesSetCutoff(option="ByLine", value=750e9)

# Use an automatic agenda
ws.propmat_clearsky_agendaSetAutomatic()

# Arts setup (No NLTE, no polarization, and standard isotopologue ratios)
ws.nlte_do = 0
ws.stokes_dim = 1
ws.isotopologue_ratiosInitFromBuiltin()

# Settings (Standard atmosphere midlatitude-summer)
ws.Touch(ws.jacobian_quantities)
ws.Touch(ws.select_abs_species)
ws.Touch(ws.rtp_nlte)
ws.rtp_mag = [10e-6, 20e-6, 40e-6]
ws.rtp_los = [45, 45]
ws.rtp_pressure = 110000
ws.rtp_temperature = 2.942000e+02
ws.f_grid = pyarts.arts.convert.kaycm2freq(np.linspace(300, 3000, 1000))
ws.rtp_vmr = [1.000869e-09, 1.000869e-14, 2.850472e-06, 1.501303e-07,
              3.019448e-08, 3.302947e-04, 1.877431e-02, 2.091960e-01, 1.701397e-06]

# Check that the calculations are OK
ws.lbl_checked = 1

# Start the explorer
ws.propmat_clearsky_agendaGUI()
