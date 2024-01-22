import pyarts
import numpy as np

# Limit the frequency range of stored absorption lines
fmin = 0
fmax = 3e12

# initialize ARTS
ws = pyarts.workspace.Workspace()

# We are setting up our absorption species to include hydrogen chloride,
# chlorine monoxide, carbon monoxide, nitrous oxide, and ozone
# To speed up calculations some, we limit ourself to absorption lines with less
# than 3 THz as their central frequency
# NOTE: Limiting the frequency range of absorption lines like this artificially
# removes absorption from the wings of stronger absorption lines that lie
# outside of the given freqeuncy range.  This migth affect the accuracy of the
# total absorption that is computed
ws.abs_speciesSet(
    species=[f"HCl-*-{fmin}-{fmax}",
             f"ClO-*-{fmin}-{fmax}",
             f"CO-*-{fmin}-{fmax}",
             f"N2O-*-{fmin}-{fmax}",
             f"O3-*-{fmin}-{fmax}"])

# Read the absorption lines.  These should be part of the arts-cata-data package
ws.abs_lines_per_speciesReadSpeciesSplitCatalog(basename="lines/")

# Use an automatic agenda
ws.propmat_clearsky_agendaAuto()

# Arts setup (No NLTE, no polarization, and standard isotopologue ratios)
ws.nlte_do = 0

# Settings (Standard atmosphere midlatitude-summer)
ws.path_point.los = [45, 45]
ws.f_grid = np.linspace(1e9, 3e12, 1000)
ws.atm_pointInit()
ws.atm_point.temperature = 2.942000e+02
ws.atm_point.pressure = 110000
ws.atm_point[ws.abs_species[0]] = 1.000869e-09
ws.atm_point[ws.abs_species[1]] = 1.000869e-14
ws.atm_point[ws.abs_species[2]] = 2.850472e-06
ws.atm_point[ws.abs_species[3]] = 1.501303e-07
ws.atm_point[ws.abs_species[4]] = 3.019448e-08
ws.atm_point.mag = [10e-6, 20e-6, 40e-6]

# Start the explorer
ws.propmat_clearsky_agendaGUI()
