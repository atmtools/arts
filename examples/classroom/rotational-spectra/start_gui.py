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
ws.absorption_speciesSet(
    species=[f"HCl-*-{fmin}-{fmax}",
             f"ClO-*-{fmin}-{fmax}",
             f"CO-*-{fmin}-{fmax}",
             f"N2O-*-{fmin}-{fmax}",
             f"O3-*-{fmin}-{fmax}"])

# Read the absorption lines.  These should be part of the arts-cata-data package
ws.AbsorptionReadSpeciesSplitCatalogs()

# Use an automatic agenda
ws.propagation_matrix_agendaAuto()

# Settings (Standard atmosphere midlatitude-summer)
ws.propagation_path_point.los = [45, 45]
ws.frequency_grid = np.linspace(1e9, 3e12, 1000)
ws.atmospheric_pointInit()
ws.atmospheric_point.temperature = 2.942000e+02
ws.atmospheric_point.pressure = 110000
ws.atmospheric_point[ws.absorption_species[0]] = 1.000869e-09
ws.atmospheric_point[ws.absorption_species[1]] = 1.000869e-14
ws.atmospheric_point[ws.absorption_species[2]] = 2.850472e-06
ws.atmospheric_point[ws.absorption_species[3]] = 1.501303e-07
ws.atmospheric_point[ws.absorption_species[4]] = 3.019448e-08
ws.atmospheric_point.mag = [10e-6, 20e-6, 40e-6]

# Start the explorer
ws.propagation_matrix_agendaGUI()
