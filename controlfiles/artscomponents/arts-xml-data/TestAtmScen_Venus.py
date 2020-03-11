#
# Testing functionality (meeting format requirements, etc.) of basic atmospheric
#  scenario data.
#
# General test setup: reading in raw data, regridding to common p-grid (1D),
#  executing standard pre-RT calc internal test method atmfields_checkedCalc.
#
#
# This case is for Venus and specifically tests
#
# (CASES A-E)
#  - the five Venus scenarions: Venus.spicav.night, Venus.spicav.night_cold,
#     Venus.vira.night, Venus.vira.day, Venus.vira.day_highlat
#  - t, z, and all abs species vmr (the ones that follow the basename convention
#     are caught by abs_speciesDefineAllInScenario and AtmRawRead; others are
#     derived manually) in the case folder AND in the "unused/" subfolder (for
#     the latter we have to tweak abs_species setting a little).
#  - regridding to a pressure grid taken from the read-in altitude grid
#
# Jana Mendrok 2013-02-26

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_venus.arts")
# 1-dimensional atmosphere
ws.AtmosphereSet1D()
ws.GriddedField3Create("gf3tmp")
ws.StringCreate("caseext")
ws.StringCreate("casefull")
ws.StringCreate("atmcase")
ws.IndexCreate("ncases")
# set basic case folder
ws.StringCreate("basename")
ws.StringSet(ws.basename, "planets/Venus/MPS/")
# Array with case names
ws.ArrayOfStringCreate("atmcasearray")
ws.ArrayOfStringSet(
    ws.atmcasearray,
    [
        "Venus.spicav.night",
        "Venus.spicav.night_cold",
        "Venus.vira.night",
        "Venus.vira.day",
        "Venus.vira.day_highlat",
    ],
)
#####
# CASES A-C (night, 10+7+2 vmr profiles) and D-E (day, 14+7+2 vmr profiles)
#####
# we go with a foorloop through the different cases
#  day & night have different number of abs_species, but identical species that
#  have to be added manually (not following standard basename or "unused/"
#  species, respectively).
@arts_agenda
def forloop_agenda(ws):
    # construct atmcase name
    ws.Copy(ws.atmcase, ws.basename)
    ws.Extract(ws.casefull, ws.atmcasearray, ws.forloop_index)
    ws.Append(ws.atmcase, ws.casefull)
    ws.StringSet(ws.caseext, "/")
    ws.Append(ws.atmcase, ws.caseext)
    ws.Append(ws.atmcase, ws.casefull)
    ws.Print(ws.atmcase)
    # derive absspecies with standard name from scenario
    ws.abs_speciesDefineAllInScenario(basename=ws.atmcase)
    # WriteXMLIndexed( "ascii", forloop_index,
    #                 abs_species, "TestAtmScen_Venus_allInScen.abs_species" )
    ws.AtmRawRead(basename=ws.atmcase)
    # adding species or variants that do not follow the general naming convention
    ws.abs_speciesAdd(species=["H2O-162"])
    ws.Copy(ws.casefull, ws.atmcase)
    ws.StringSet(ws.caseext, ".H2O-162_high")
    ws.Append(ws.casefull, ws.caseext)
    ws.ReadXML(ws.gf3tmp, ws.casefull)
    ws.Append(ws.vmr_field_raw, ws.gf3tmp)
    ws.abs_speciesAdd(species=["H2O-162"])
    ws.Copy(ws.casefull, ws.atmcase)
    ws.StringSet(ws.caseext, ".H2O-162_low")
    ws.Append(ws.casefull, ws.caseext)
    ws.ReadXML(ws.gf3tmp, ws.casefull)
    ws.Append(ws.vmr_field_raw, ws.gf3tmp)
    ws.abs_speciesAdd(species=["H2O-162"])
    ws.Copy(ws.casefull, ws.atmcase)
    ws.StringSet(ws.caseext, ".H2O-162_mid")
    ws.Append(ws.casefull, ws.caseext)
    ws.ReadXML(ws.gf3tmp, ws.casefull)
    ws.Append(ws.vmr_field_raw, ws.gf3tmp)
    ws.abs_speciesAdd(species=["H2O-162"])
    ws.Copy(ws.casefull, ws.atmcase)
    ws.StringSet(ws.caseext, ".H2O-162_uncorrected")
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
    ws.abs_speciesAdd(species=["H2O"])
    ws.Copy(ws.casefull, ws.atmcase)
    ws.StringSet(ws.caseext, ".H2O_mid")
    ws.Append(ws.casefull, ws.caseext)
    ws.ReadXML(ws.gf3tmp, ws.casefull)
    ws.Append(ws.vmr_field_raw, ws.gf3tmp)
    # reset atmcase name for profiles in "unused/" folder
    ws.Copy(ws.atmcase, ws.basename)
    ws.Extract(ws.casefull, ws.atmcasearray, ws.forloop_index)
    ws.Append(ws.atmcase, ws.casefull)
    ws.StringSet(ws.caseext, "/unused/")
    ws.Append(ws.atmcase, ws.caseext)
    ws.Append(ws.atmcase, ws.casefull)
    ws.Print(ws.atmcase)
    # add the profiles (and a dummy abs_species) of the unused species
    ws.abs_speciesAdd(species=["N2"])
    ws.Copy(ws.casefull, ws.atmcase)
    ws.StringSet(ws.caseext, ".S2")
    ws.Append(ws.casefull, ws.caseext)
    ws.ReadXML(ws.gf3tmp, ws.casefull)
    ws.Append(ws.vmr_field_raw, ws.gf3tmp)
    ws.abs_speciesAdd(species=["N2"])
    ws.Copy(ws.casefull, ws.atmcase)
    ws.StringSet(ws.caseext, ".SO3")
    ws.Append(ws.casefull, ws.caseext)
    ws.ReadXML(ws.gf3tmp, ws.casefull)
    ws.Append(ws.vmr_field_raw, ws.gf3tmp)
    # now derive common p_grid and regrid atm fields to this
    ws.p_gridFromZRaw(ws.p_grid, ws.z_field_raw, 0)
    ws.AtmFieldsCalc(vmr_zeropadding=1)
    ws.atmfields_checkedCalc(bad_partition_functions_ok=1)
    # WriteXML( "ascii", p_grid )
    # WriteXML( "ascii", z_field )
    # WriteXML( "ascii", t_field )
    # WriteXML( "ascii", vmr_field_raw )
    # WriteXMLIndexed( "ascii", forloop_index, vmr_field )


ws.forloop_agenda = forloop_agenda

ws.nelemGet(ws.ncases, ws.atmcasearray)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
