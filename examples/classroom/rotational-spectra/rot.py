import pyarts
import numpy as np

# Limit the frequency range of stored absorption lines
fmin = 0
fmax = 3e12

ws = pyarts.workspace.Workspace()

# Species setup and catalog reading
ws.abs_speciesSet(
    species=[f"HCl-*-{fmin}-{fmax}",
             f"ClO-*-{fmin}-{fmax}",
             f"CO-*-{fmin}-{fmax}",
             f"N2O-*-{fmin}-{fmax}",
             f"O3-*-{fmin}-{fmax}"])
ws.abs_lines_per_speciesReadSpeciesSplitCatalog(basename="lines/")

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
ws.f_grid = np.linspace(1e9, 3e12, 1000)
ws.rtp_vmr = [1.000869e-09, 1.000869e-14,
              2.850472e-06, 1.501303e-07, 3.019448e-08]

# Check that the calculations are OK
ws.lbl_checkedCalc()

# Start the explorer
ws.propmat_clearsky_agendaGUI()
