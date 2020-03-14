# DEFINITIONS:  -*-sh-*-
# This is a test doing simulations for MSG-ICI instrument.
#
# So far, it's doing clearsky-pencilbeam simulations from and for (planned) ICI
#  orbit and observation geometry considering planned radiometer band
#  characteristics.
#
# Author: Jana Mendrok

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
# Agenda for scalar gas absorption calculation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
# Basic settings (already needed in sensor part)
# ---
ws.AtmosphereSet1D()
ws.IndexSet(ws.stokes_dim, 1)
# 1. General
# -----------
ws.output_file_formatSetZippedAscii()
# AMSU uses Planck brightness temperatures
# ---
ws.StringSet(ws.iy_unit, "PlanckBT")
# AMSU uses Planck brightness temperatures
# ---
ws.StringSet(ws.iy_unit, "PlanckBT")
# (standard) emission calculation
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
# cosmic background radiation
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
# standard surface agenda (i.e., make use of surface_rtprop_agenda)
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
# sensor-only path
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
# no refraction
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
#
# modifiy the maximum propagation step, from the default to 250 m :
#
ws.NumericSet(ws.ppath_lmax, 250.0)
# Surface
ws.Copy(
    ws.surface_rtprop_agenda,
    ws.surface_rtprop_agenda__Specular_NoPol_ReflFix_SurfTFromt_surface,
)
# 2. Spectroscopy
# ----------------
# We take a smaller cutoff, since the line-by-line calculation is
# only for O3, where only the local lines matter.
# Could be speed-optimized further by selecting only the relevant
# lines from the line list.
ws.abs_speciesSet(species=["H2O-PWR98", "O2-PWR93", "N2-SelfContStandardType"])
# Read HITRAN catalog (needed for O3):
# abs_linesReadFromHitran( abs_lines,
#                         "/storage3/data/catalogue/hitran/hitran2012/HITRAN2012.par",
#                         150e9,
#                         700e9 )
# WriteXML( "ascii", abs_lines, "ici.hitran12_lines.xml" )
# ReadXML( abs_lines, "ici.hitran12_lines.xml" )
# abs_lines_per_speciesCreateFromLines
ws.abs_lines_per_speciesSetEmpty()
# WARNING: If you redefine abs_species, and want to do a line-by-line
# calculation, you also have to call
# abs_lines_per_speciesCreateFromLines again.
# 3. Sensor:
# -----------
ws.execute_controlfile("instruments/ici/ici_fast.arts")
# 4. Atmosphere
# --------------
# Atmospheric profiles are stored in an ArrayOfGriddedField4.
# It contains one GriddedField4 for each atmospheric state.
#
ws.ReadXML(ws.batch_atm_fields_compact, "../../testdata/chevallierl91_all_extract.xml")
# add constant profiles for O2 and N2
ws.batch_atm_fields_compactAddConstant(name="abs_species-O2", value=0.2095)
ws.batch_atm_fields_compactAddConstant(name="abs_species-N2", value=0.7808)
# Set parameters for lookup table
# ---
# Arguments omitted for better maintainability of this test file.
# abs_lookupSetupWide
ws.abs_lookupSetupBatch()
# Create the lookup table
# ---
ws.abs_xsec_agenda_checkedCalc()
ws.lbl_checkedCalc()
ws.jacobianOff()
ws.abs_lookupCalc()
# Set propmat_clearsky_agenda to use lookup table
# ---
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__LookUpTable)
# Set up RT calculation
# =====================
# Set surface reflectivity
# ---
# Here we take a value representative for the sea surface.
# NumericSet( surface_emissivity, 0.6 )  <--- Old, replaced by:
ws.VectorSetConstant(ws.surface_scalar_reflectivity, 1, 0.4)
# Set the agenda for batch calculations:
# ---
#
@arts_agenda
def ybatch_calc_agenda(ws):
    # Extract the atmospheric profiles for this case:
    ws.Extract(ws.atm_fields_compact, ws.batch_atm_fields_compact, ws.ybatch_index)
    # Split up *atm_fields_compact* to
    # generate p_grid, t_field, z_field, vmr_field:
    ws.AtmFieldsAndParticleBulkPropFieldFromCompact()
    # Optionally set Jacobian parameters.
    # uncomment this for NO jacobian calculations
    ws.jacobianOff()
    # Uncomment this block if you want Jacobians. Attention, it slows down the
    # computation a lot.
    # Also, you can add other Jacobians here, for example for temperature.
    #  jacobianInit
    #  jacobianAddAbsSpecies( jacobian_quantities, jacobian_agenda,
    #                         atmosphere_dim,
    #                         p_grid, lat_grid, lon_grid,
    #                         p_grid, lat_grid, lon_grid,
    #                         "H2O-PWR98",
    #                         "rel")
    #  jacobianClose
    # No scattering
    ws.cloudboxOff()
    # get some surface properties from corresponding atmospheric fields
    ws.Extract(ws.z_surface, ws.z_field, 0)
    ws.Extract(ws.t_surface, ws.t_field, 0)
    # Perform RT calculations
    # ---
    ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
    ws.atmgeom_checkedCalc()
    ws.cloudbox_checkedCalc()
    ws.sensor_checkedCalc()
    ws.yCalc()


ws.ybatch_calc_agenda = ybatch_calc_agenda

# Set number of batch cases:
ws.nelemGet(ws.ybatch_n, ws.batch_atm_fields_compact)
# IndexSet(ybatch_start, 2)
# IndexSet(ybatch_n, 2)
# Execute the batch calculations:
# ---
ws.propmat_clearsky_agenda_checkedCalc()
ws.ybatchCalc()
# Store result matrix:
# ---
ws.WriteXML("ascii", ws.ybatch)
# WriteXML( "ascii", ybatch_jacobians )
ws.ArrayOfVectorCreate("ybatch_ref")
ws.ReadXML(ws.ybatch_ref, "TestICI_fast.ybatch.ref.xml")
ws.Compare(
    ws.ybatch,
    ws.ybatch_ref,
    0.2,
    "Total radiance should be close to the reference values",
)
