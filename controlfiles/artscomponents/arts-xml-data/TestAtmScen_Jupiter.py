#
# Testing functionality (meeting format requirements, etc.) of basic atmospheric
#  scenario data.
#
# General test setup: reading in raw data, regridding to common p-grid (1D),
#  executing standard pre-RT calc internal test method atmfields_checkedCalc.
#
#
# This case is for Jupiter and specifically tests
#
# (CASE A)
#  - from MPS/Jupiter.mean: t, z, and all abs species vmr (the ones that follow
#     the basename convention are caught by abs_speciesDefineAllInScenario and
#     AtmRawRead; others are derived manually) in the case folder AND in the
#     "unused/" subfolder (for the latter we have to tweak abs_species setting a
#     little).
# (CASE B)
# - from MPS/Jupiter.oval: t, z (re-using the vmr profiles from Jupiter.mean).
#
# - regridding to a pressure grid taken from the read-in altitude grid
# - regridding to a user defined pressure grid (we chose one that is covers a
#    slightly larger p-range than the original data in order to test the
#    extrapolation and zero-padding features).
#
# Jana Mendrok 2013-02-26

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_jupiter.arts")
# 1-dimensional atmosphere
ws.AtmosphereSet1D()
ws.GriddedField3Create("gf3tmp")
ws.StringCreate("caseext")
ws.StringCreate("casefull")
ws.StringCreate("atmcase")
#####
# CASE A
#####
# set atmospheric scenario
ws.StringSet(ws.atmcase, "planets/Jupiter/MPS/Jupiter.mean/Jupiter.mean")
# derive abs species from scenario data
ws.abs_speciesDefineAllInScenario(basename=ws.atmcase)
# WriteXML( "ascii", abs_species, "TestAtmScen_Jupiter_allInScen.abs_species.xml" )
# get atm scenario raw data
ws.AtmRawRead(basename=ws.atmcase)
# adding species or variants that do not follow the general naming convention
ws.abs_speciesAdd(species=["CH4-212"])
ws.Copy(ws.casefull, ws.atmcase)
ws.StringSet(ws.caseext, ".CH4-212")
ws.Append(ws.casefull, ws.caseext)
ws.ReadXML(ws.gf3tmp, ws.casefull)
ws.Append(ws.vmr_field_raw, ws.gf3tmp)
ws.abs_speciesAdd(species=["H2-12"])
ws.Copy(ws.casefull, ws.atmcase)
ws.StringSet(ws.caseext, ".H2-12")
ws.Append(ws.casefull, ws.caseext)
ws.ReadXML(ws.gf3tmp, ws.casefull)
ws.Append(ws.vmr_field_raw, ws.gf3tmp)
ws.abs_speciesAdd(species=["H2O"])
ws.Copy(ws.casefull, ws.atmcase)
ws.StringSet(ws.caseext, ".H2O_high")
ws.Append(ws.casefull, ws.caseext)
ws.ReadXML(ws.gf3tmp, ws.casefull)
ws.Append(ws.vmr_field_raw, ws.gf3tmp)
ws.abs_speciesAdd(species=["H2O"])
ws.Copy(ws.casefull, ws.atmcase)
ws.StringSet(ws.caseext, ".H2O_low")
ws.Append(ws.casefull, ws.caseext)
ws.ReadXML(ws.gf3tmp, ws.casefull)
ws.Append(ws.vmr_field_raw, ws.gf3tmp)
ws.abs_speciesAdd(species=["HCN"])
ws.Copy(ws.casefull, ws.atmcase)
ws.StringSet(ws.caseext, ".HCN_upperlim")
ws.Append(ws.casefull, ws.caseext)
ws.ReadXML(ws.gf3tmp, ws.casefull)
ws.Append(ws.vmr_field_raw, ws.gf3tmp)
ws.abs_speciesAdd(species=["NH3"])
ws.Copy(ws.casefull, ws.atmcase)
ws.StringSet(ws.caseext, ".NH3_high")
ws.Append(ws.casefull, ws.caseext)
ws.ReadXML(ws.gf3tmp, ws.casefull)
ws.Append(ws.vmr_field_raw, ws.gf3tmp)
ws.abs_speciesAdd(species=["NH3"])
ws.Copy(ws.casefull, ws.atmcase)
ws.StringSet(ws.caseext, ".NH3_low")
ws.Append(ws.casefull, ws.caseext)
ws.ReadXML(ws.gf3tmp, ws.casefull)
ws.Append(ws.vmr_field_raw, ws.gf3tmp)
# WriteXML( "ascii", abs_species, "TestAtmScen_Jupiter_allvalid.abs_species.xml" )
#####
# we also test the unused species for correct format. as these species are not
# valid ARTS species, we have to assign those profile data to some other abs
# species. we take N2 (but any would be ok, as long as we don't do any abs calc
# with it.
# reset atmospheric scenario
ws.StringSet(ws.atmcase, "planets/Jupiter/MPS/Jupiter.mean/unused/Jupiter.mean")
ws.abs_speciesAdd(species=["N2"])
ws.Copy(ws.casefull, ws.atmcase)
ws.StringSet(ws.caseext, ".AsH3")
ws.Append(ws.casefull, ws.caseext)
ws.ReadXML(ws.gf3tmp, ws.casefull)
ws.Append(ws.vmr_field_raw, ws.gf3tmp)
ws.abs_speciesAdd(species=["N2"])
ws.Copy(ws.casefull, ws.atmcase)
ws.StringSet(ws.caseext, ".C6H6")
ws.Append(ws.casefull, ws.caseext)
ws.ReadXML(ws.gf3tmp, ws.casefull)
ws.Append(ws.vmr_field_raw, ws.gf3tmp)
ws.abs_speciesAdd(species=["N2"])
ws.Copy(ws.casefull, ws.atmcase)
ws.StringSet(ws.caseext, ".GeH4")
ws.Append(ws.casefull, ws.caseext)
ws.ReadXML(ws.gf3tmp, ws.casefull)
ws.Append(ws.vmr_field_raw, ws.gf3tmp)
ws.abs_speciesAdd(species=["N2"])
ws.Copy(ws.casefull, ws.atmcase)
ws.StringSet(ws.caseext, ".H")
ws.Append(ws.casefull, ws.caseext)
ws.ReadXML(ws.gf3tmp, ws.casefull)
ws.Append(ws.vmr_field_raw, ws.gf3tmp)
#####
# A-1) p_grid initialized from given altitude grid
#####
# we need to include negative altitudes!
ws.p_gridFromZRaw(ws.p_grid, ws.z_field_raw, 0)
ws.AtmFieldsCalc(vmr_zeropadding=1)
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
ws.atmgeom_checkedCalc()
# WriteXML( "ascii", p_grid )
# WriteXML( "ascii", z_field )
# WriteXML( "ascii", t_field )
# WriteXML( "ascii", vmr_field_raw )
# WriteXML( "ascii", vmr_field )
#####
# A-2) p_grid set to a user defined grid (surely requries interpolation to calc-grid(s))
#####
ws.VectorNLogSpace(ws.p_grid, 401, 1020000.0, 9.8e-10)
ws.AtmFieldsCalc(vmr_zeropadding=1)
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
ws.atmgeom_checkedCalc()
# WriteXML( "ascii", p_grid )
# WriteXML( "ascii", z_field )
# WriteXML( "ascii", t_field )
# WriteXML( "ascii", vmr_field_raw )
# WriteXML( "ascii", vmr_field )
#####
# CASE B
#####
# Now doing the second case folder. however, this does not include any abs
# species (or vmr) data, but it does have t and z data, which we like to test.
# Hence, we only replace those data and keep the raw vmr profiles from before.
# As we don't read a full scenario (and don't want to overwrite the previous vrm
# raw data, we ReadXML the t and z data.
# set atmospheric scenario
ws.StringSet(ws.atmcase, "planets/Jupiter/MPS/Jupiter.oval/Jupiter.oval")
ws.Copy(ws.casefull, ws.atmcase)
ws.StringSet(ws.caseext, ".t")
ws.Append(ws.casefull, ws.caseext)
ws.ReadXML(ws.t_field_raw, ws.casefull)
ws.Copy(ws.casefull, ws.atmcase)
ws.StringSet(ws.caseext, ".z")
ws.Append(ws.casefull, ws.caseext)
ws.ReadXML(ws.z_field_raw, ws.casefull)
#####
# B-1) p_grid initialized from given altitude grid
#####
# we need to include negative altitudes!
ws.p_gridFromZRaw(ws.p_grid, ws.z_field_raw, 0)
ws.AtmFieldsCalc(vmr_zeropadding=1)
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
ws.atmgeom_checkedCalc()
# WriteXML( "ascii", p_grid )
# WriteXML( "ascii", z_field )
# WriteXML( "ascii", t_field )
# WriteXML( "ascii", vmr_field_raw )
# WriteXML( "ascii", vmr_field )
