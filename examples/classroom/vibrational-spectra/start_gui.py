import pyarts
import numpy as np

# initialize ARTS
ws = pyarts.workspace.Workspace()

# These are the same species as for the rotational spectra (see
# ../rotational-spectra/start_gui.py) but we also add water, carbon dioxide
# and methane since these are sometimes more interesting in the infrared.
# Gone are also the frequency limits, we will later use another speed-up
# technique to achieve faster copmpute times
ws.abs_speciesSet(
    species=["HCl",
             "ClO",
             "CO",
             "N2O",
             "O3",
             "H2O",
             "CO2",
             "CH4"])

# Read the absorption lines.  These should be part of the arts-cata-data package
ws.abs_lines_per_speciesReadSpeciesSplitCatalog(basename="lines/")

# Speed-up technique of using a cutoff frequency
# This will apply a cutoff of each absorption line at a frequency of 750 GHz
# away from the pressure shifted line center.  The value at the cutoff frequency
# is subsequently used to remove the value from within the line center.
# the latter modification is in line with how water continua models deal with
# some of its unknown or just difficult to model absorption features, however
# we will not consider this absorption in the current example
# NOTE: Using a cutoff frequency across all species like this degrades the
# quality of the absorption calculations, which can be seen by changing the
# Y-scale of the absorption to logarithmic
for bands in ws.abs_lines_per_species:
    for band in bands:
        band.cutoff = pyarts.arts.options.AbsorptionCutoffType.ByLine
        band.cutofffreq = 750e9

# Use an automatic agenda
ws.propagation_matrix_agendaAuto()

# Arts setup (No NLTE, no polarization, and standard isotopologue ratios)
ws.nlte_do = 0

# Settings (Standard atmosphere midlatitude-summer)
ws.path_point.los = [45, 45]
ws.frequency_grid = pyarts.arts.convert.kaycm2freq(np.linspace(300, 3000, 1000))
ws.atmospheric_pointInit()
ws.atmospheric_point.temperature = 2.942000e+02
ws.atmospheric_point.pressure = 110000
ws.atmospheric_point[ws.abs_species[0]] = 1.000869e-09
ws.atmospheric_point[ws.abs_species[1]] = 1.000869e-14
ws.atmospheric_point[ws.abs_species[2]] = 2.850472e-06
ws.atmospheric_point[ws.abs_species[3]] = 1.501303e-07
ws.atmospheric_point[ws.abs_species[4]] = 3.019448e-08
ws.atmospheric_point[ws.abs_species[5]] = 1.877431e-02
ws.atmospheric_point[ws.abs_species[6]] = 3.302947e-04
ws.atmospheric_point[ws.abs_species[7]] = 1.701397e-06
ws.atmospheric_point.mag = [10e-6, 20e-6, 40e-6]

# Start the explorer
ws.propagation_matrix_agendaGUI()
