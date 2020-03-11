#
# Demonstration and test of adjusting z_field according to hydrostatic equilibrium
#
# Jana Mendrok 2013-01-29

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
# Number of Stokes components to be computed
#
ws.IndexSet(ws.stokes_dim, 4)
# Frequency grid
#
ws.VectorNLogSpace(ws.f_grid, 101, 100000000.0, 5000000000.0)
# Dimensionality of the atmosphere
#
ws.AtmosphereSet1D()
#
# Some settings resulting in no gas absorption
#
# on-the-fly absorption
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
# Definition of species
#
ws.abs_speciesSet(species=["H2O", "N2", "O2"])
# No line data needed here
#
ws.abs_lines_per_speciesSetEmpty()
# A pressure grid rougly matching 0 to 1000 km, in steps of 2km.
#
ws.VectorNLogSpace(ws.p_grid, 501, 101300.0, 1e-80)
# Surface altitude
ws.MatrixSetConstant(ws.z_surface, 1, 1, 0.0)
# Neutral atmosphere, up to about 95 km
ws.AtmRawRead(basename="testdata/tropical")
# Tempature and z also covering the ionosphere
#
ws.ReadXML(ws.t_field_raw, "testdata/tropical.expanded.t.xml")
ws.ReadXML(ws.z_field_raw, "testdata/tropical.expanded.z.xml")
# Interpolate to p_grid (VMR is "zero padded")
ws.AtmFieldsCalc(vmr_zeropadding=1)
# Apply HSE
#
# (this roughly recreates the original altitudes for the input electron density
#  and magnetic field)
#
ws.VectorSet(ws.lat_true, array([0.0]))
ws.VectorSet(ws.lon_true, array([0.0]))
ws.NumericSet(ws.p_hse, 101300.0)
ws.NumericSet(ws.z_hse_accuracy, 10.0)
#
ws.atmfields_checkedCalc()
ws.z_fieldFromHSE()
# No jacobian calculations
#
ws.jacobianOff()
# No scattering
#
ws.cloudboxOff()
# Check model atmosphere
#
ws.atmfields_checkedCalc()
ws.cloudbox_checkedCalc()
# WriteXML( "ascii", z_field, "z_fieldFromHSE_REFERENCE.xml" )
# Expected results
#
ws.Tensor3Create("zREFERENCE")
#
ws.ReadXML(ws.zREFERENCE, "z_fieldFromHSE_REFERENCE.xml")
# Check
#
ws.Compare(ws.z_field, ws.zREFERENCE, 0.0001)
# SIDE ISSUE: testing effect of H2O for z_fieldFromHSE
# (result: dev<50m for z<1400km - basically indicates that H2O is not that important...)
#
# we need to re-prepare the atmo setup without H2O
# abs_speciesSet( species=
#                ["N2",
#                 "O2"] )
# AtmRawRead( t_field_raw, z_field_raw, vmr_field_raw, abs_species,
#            "testdata/tropical" )
# ReadXML( t_field_raw, "testdata/tropical.expanded.t.xml" )
# ReadXML( z_field_raw, "testdata/tropical.expanded.z.xml" )
# AtmFieldsCalc( t_field, z_field, vmr_field, p_grid,
#               lat_grid, lon_grid, t_field_raw, z_field_raw,
#               vmr_field_raw, atmosphere_dim )
# atmfields_checkedCalc
# atmgeom_checkedCalc
# z_fieldFromHSE
# WriteXML( "ascii", z_field, "z_fieldFromHSE_noH2O.xml" )
