# DEFINITIONS:  -*-sh-*-

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general")
ws.execute_controlfile("general/continua")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__LookUpTable)
# Create some variables:
ws.NumericCreate("fmin")
ws.NumericCreate("fmax")
ws.NumericCreate("wvl_min")
ws.NumericCreate("wvl_max")
ws.VectorCreate("wvl_grid")
# Set minimum and maximum wavelength
ws.NumericSet(ws.wvl_min, 1e-05)
# 10.00 um
ws.NumericSet(ws.wvl_max, 1.001e-05)
# 10.01 um
# Convert to Hz, maximum wavelength = minimum frequency:
ws.FrequencyFromWavelength(ws.fmax, ws.wvl_min)
ws.FrequencyFromWavelength(ws.fmin, ws.wvl_max)
ws.Print(ws.fmin, 1)
ws.Print(ws.fmax, 1)
# Define equidistant frequency grid (1000 grid points):
ws.VectorNLinSpace(ws.f_grid, 1000, ws.fmin, ws.fmax)
# Read HITRAN data (for this example we use a reduced test catalog)
ws.abs_linesReadFromArts(ws.abs_lines, "testdata/abs_lines_IR.xml.gz", ws.fmin, ws.fmax)
# Set species to be considered in line-by-line calculation
ws.abs_speciesSet(
    species=[
        "H2O, H2O-SelfContCKDMT100, H2O-ForeignContCKDMT100",
        "CO2, CO2-CKDMT100",
        "O3",
        "N2O",
        "O2, O2-CIAfunCKDMT100",
        "HNO3",
        "N2, N2-CIAfunCKDMT100, N2-CIArotCKDMT100",
    ]
)
# Alternatively select all species that we can find in the scenario:
# abs_speciesDefineAllInScenario( basename="testdata/tropical" )
# This separates the lines into the different tag groups and creates
# the workspace variable `abs_lines_per_species':
ws.abs_lines_per_speciesCreateFromLines()
# Dimensionality of the atmosphere
#
ws.AtmosphereSet1D()
# Atmospheric profiles (there is several data in the
# arts-xml-data/atmosphere directory, fascod includes the standard
# atmospheres which we also have in libRadtran (altitude only up to
# 95 km !!). When you want to use the molecular_tau_file from arts,
# the atmosphere_file for uvspec must correspond to the ARTS
# atmosphere files which are defined here!!)
ws.AtmRawRead(basename="testdata/tropical")
# Extract pressure grid from atmosphere files (this is the vertical
# coordinate for all calculations, can be specified as you like)
ws.p_gridFromZRaw(ws.p_grid, ws.z_field_raw)
ws.VectorSet(ws.lat_grid, array([], dtype=float64))
ws.VectorSet(ws.lon_grid, array([], dtype=float64))
# Now interpolate all the raw atmospheric input onto the pressure
# grid and create the atmospheric variables `t_field', `z_field', `vmr_field'
ws.AtmFieldsCalc()
# Initialize the input variables of abs_coefCalc from the Atm fields:
ws.AbsInputFromAtmFields()
# Non-linear species
ws.abs_speciesSet(abs_species=ws.abs_nls, species=[])
# Perturbation if lookup-table should be created that can be used for a wide range of atmospheric conditions
ws.VectorSet(ws.abs_t_pert, array([], dtype=float64))
ws.VectorSet(ws.abs_nls_pert, array([], dtype=float64))
ws.abs_xsec_agenda_checkedCalc()
ws.jacobianOff()
# If you want to do calculations for various atmopheric conditions, check out how the lookup table works.
# This saves a lot of time!!!
ws.abs_lookupCalc()
# output_file_formatSetBinary
# WriteXML (output_file_format, abs_lookup )
# Calculate absorption field:
ws.IndexSet(ws.f_index, -1)
# calculate all frequencies
ws.IndexSet(ws.stokes_dim, 1)
ws.atmfields_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.propmat_clearsky_fieldCalc()
# Write molecular_tau_file for libRadtran
ws.WriteMolTau(
    ws.f_grid, ws.z_field, ws.propmat_clearsky_field, ws.atmosphere_dim, "TestMolTau.nc"
)
