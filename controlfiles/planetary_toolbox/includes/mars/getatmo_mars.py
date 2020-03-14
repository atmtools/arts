################################################################################
#                                                                              #
# Unless further variables or options for existing variables are introduced,   #
# DO NOT MODIFY this file! This is only a helper file!                         #
#                                                                              #
################################################################################
#                                                                              #
# This file does the actual work of selecting and reading in the RAW           #
# atmosphere data for Mars as specified by the user. For user specification    #
# use, e.g., DemoMarsAtmo1D.arts (or its 3D equivalent) as template. The       #
# template also contains the detailed information on which indices are linked  #
# to which specific value/selection for each of the variables. The full        #
# arrays, which the indices refer to and from which the actual values are      #
# extracted, are defined in atmo_mars.arts (hence, atmo_mars.arts needs to be  #
# included before the file at hand).                                           #
#                                                                              #
# This file expects the following input parameters:                            #
#   Ls             (Index)           The season of the atmospheric scenario.   #
#   daytime        (Index)           The daytime of the atmospheric scenario.  #
#   dust           (Index)           The dust loading of the atmospheric       #
#                                     scenario.                                #
#   solar          (Index)           The solar activity level of the           #
#                                     atmospheric scenario.                    #
#   basespecies    (ArrayOfIndex)    The abs_species to use (includes only     #
#                                     such with on/off options only).          #
#   h2ospecies     (ArrayOfIndex)    H2O setup selected (off/low/medium/high). #
#   hdospecies     (ArrayOfIndex)    HDO setup selected                        #
#                                     (off/low/medium/high/uncorrected).       #
#   Necase         (ArrayOfIndex)    Electron density setup selected           #
#                                     (off/different sza-dependencies).        #
#                                                                              #
# Files to be included before this file:                                       #
#   includes/mars/atmo_mars.arts                                               #
#   includes/common/createvars.arts                                            #
#                                                                              #
# It provides following output:                                                #
#   z_field_raw        as the WSV                                              #
#   t_field_raw        as the WSV                                              #
#   vmr_field_raw      as the WSV                                              #
#   abs_species        as the WSV                                              #
#                                                                              #
################################################################################

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# We will need to dummy-store some data in files to be able to export data from
# forloops. So we create some dummy names.
ws.StringSet(ws.tmpformat, "binary")
ws.StringSet(ws.vmrtmp, "vmrtmp_mars.xml")
ws.StringSet(ws.abstmp, "abstmp_mars.xml")
# Get data for abs_species, where filenames do follow the basic convention
# (filename=speciesname.xml)
ws.AgendaCreate("speciesloop_agenda")


@arts_agenda
def speciesloop_agenda(ws):
    ws.ReadXML(out=ws.vmr_field_raw, filename=ws.vmrtmp)
    ws.Extract(ws.strtmp, ws.casearray, ws.forloop_index)
    ws.Append(ws.specfilename, ws.strtmp)
    #  Print( specfilename, 0 )
    ws.ReadXML(ws.gf3tmp, ws.specfilename)
    ws.Append(ws.vmr_field_raw, ws.gf3tmp)
    ws.WriteXML(ws.tmpformat, ws.vmr_field_raw, ws.vmrtmp, 0)


ws.speciesloop_agenda = speciesloop_agenda

# Get data for abs_species, where filenames do NOT follow the basic convention
ws.AgendaCreate("subspeciesloop_agenda")


@arts_agenda
def subspeciesloop_agenda(ws):
    ws.ReadXML(out=ws.vmr_field_raw, filename=ws.vmrtmp)
    ws.ReadXML(out=ws.abs_species, filename=ws.abstmp)
    ws.abs_speciesAdd(species=ws.speciesname)
    ws.Extract(ws.strtmp, ws.casearray, ws.forloop_index)
    ws.Append(ws.specfilename, ws.strtmp)
    #  Print( specfilename, 0 )
    ws.ReadXML(ws.gf3tmp, ws.specfilename)
    ws.Append(ws.vmr_field_raw, ws.gf3tmp)
    ws.WriteXML(ws.tmpformat, ws.vmr_field_raw, ws.vmrtmp, 0)
    ws.WriteXML(ws.tmpformat, ws.abs_species, ws.abstmp, 0)


ws.subspeciesloop_agenda = subspeciesloop_agenda

# Read the atmospheric raw data
# ---
# first, create the casename string down to the common filename part in the
# scenario folder. Data is located in:
# Mars.Ls.daytime.dust/
#   Mars.Ls.daytime.dust.solar/
#     Mars.Ls.daytime.dust.solar.[variable.xml]
ws.Copy(ws.superatmo, ws.atmobase)
# construct upper level path name (Mars.Ls.daytime.dust)
ws.StringSet(ws.atmostr, "Mars")
ws.StringSet(ws.strtmp, ".")
ws.Append(ws.atmostr, ws.strtmp)
ws.Extract(ws.subatmo, ws.Lsarray, ws.Ls)
ws.Append(ws.atmostr, ws.subatmo)
ws.Append(ws.atmostr, ws.strtmp)
ws.Extract(ws.subatmo, ws.daytimearray, ws.daytime)
ws.Append(ws.atmostr, ws.subatmo)
ws.Append(ws.atmostr, ws.strtmp)
ws.Extract(ws.subatmo, ws.dustarray, ws.dust)
ws.Append(ws.atmostr, ws.subatmo)
# append upper level path name (Mars.Ls.daytime.dust) to base path
ws.Append(ws.superatmo, ws.atmostr)
ws.StringSet(ws.strtmp, "/")
ws.Append(ws.superatmo, ws.strtmp)
# go on to construct lower level path name (Mars.Ls.daytime.dust.solar)
ws.StringSet(ws.strtmp, ".")
ws.Append(ws.atmostr, ws.strtmp)
ws.Extract(ws.subatmo, ws.solararray, ws.solar)
ws.Append(ws.atmostr, ws.subatmo)
# append lower level path name (Mars.Ls.daytime.dust.solar) to base/upper-level
# path
ws.Append(ws.superatmo, ws.atmostr)
ws.StringSet(ws.strtmp, "/")
ws.Append(ws.superatmo, ws.strtmp)
# append base file name (Mars.Ls.daytime.dust.solar.) to path construction
ws.Append(ws.superatmo, ws.atmostr)
ws.StringSet(ws.strtmp, ".")
ws.Append(ws.superatmo, ws.strtmp)
# eventually, copy full path-basefilename combi to atmostr, to be consistent
# with setups for the other planets, where we used atmostr for the atmo scenario
# file location
ws.Copy(ws.atmostr, ws.superatmo)
ws.StringSet(ws.infostr, "Atmospheric data taken from: ")
ws.Append(ws.infostr, ws.atmostr)
ws.Print(ws.infostr)
# second, we construct the name for the specific data files one-by-one and read
# into corresponding variable
ws.Touch(ws.vmr_field_raw)
ws.Touch(ws.abs_species)
ws.WriteXML(ws.tmpformat, ws.vmr_field_raw, ws.vmrtmp, 0)
ws.WriteXML(ws.tmpformat, ws.abs_species, ws.abstmp, 0)
# (1) z = Altitude
ws.Copy(ws.specfilename, ws.atmostr)
ws.StringSet(ws.strtmp, "z.xml")
ws.Append(ws.specfilename, ws.strtmp)
ws.ReadXML(ws.z_field_raw, ws.specfilename)
# (2) t = Temperature
ws.Copy(ws.specfilename, ws.atmostr)
ws.StringSet(ws.strtmp, "t.xml")
ws.Append(ws.specfilename, ws.strtmp)
ws.ReadXML(ws.t_field_raw, ws.specfilename)
# (3) Ne = electron density
ws.Copy(ws.specfilename, ws.atmostr)
ws.ArrayOfStringSet(ws.speciesname, ["free_electrons"])
ws.Select(ws.casearray, ws.Nearray, ws.Necase)
ws.nelemGet(ws.ncases, ws.casearray)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.Copy(ws.forloop_agenda, ws.subspeciesloop_agenda)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
ws.ReadXML(out=ws.abs_species, filename=ws.abstmp)
# (4) base-vmr (species without subscenarios)
ws.Copy(ws.specfilename, ws.atmostr)
ws.Select(ws.speciesname, ws.basespeciesarray, ws.basespecies)
ws.abs_speciesAdd(species=ws.speciesname)
ws.WriteXML(ws.tmpformat, ws.abs_species, ws.abstmp, 0)
ws.Select(ws.casearray, ws.basespeciesnamesarray, ws.basespecies)
ws.nelemGet(ws.ncases, ws.casearray)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.Copy(ws.forloop_agenda, ws.speciesloop_agenda)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
# (5) CH4
ws.Copy(ws.specfilename, ws.atmostr)
ws.ArrayOfStringSet(ws.speciesname, ["CH4"])
ws.Select(ws.casearray, ws.CH4array, ws.ch4species)
ws.nelemGet(ws.ncases, ws.casearray)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.Copy(ws.forloop_agenda, ws.subspeciesloop_agenda)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
ws.ReadXML(out=ws.abs_species, filename=ws.abstmp)
# (6) H2O
ws.Copy(ws.specfilename, ws.atmostr)
ws.Select(ws.speciesname, ws.H2Oarray, ws.h2ospecies)
ws.abs_speciesAdd(species=ws.speciesname)
ws.WriteXML(ws.tmpformat, ws.abs_species, ws.abstmp, 0)
ws.Select(ws.casearray, ws.H2Onamesarray, ws.h2ospecies)
ws.nelemGet(ws.ncases, ws.casearray)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.Copy(ws.forloop_agenda, ws.speciesloop_agenda)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
# now we're ready with the abs_species (and vmr_fields).
ws.ReadXML(out=ws.abs_species, filename=ws.abstmp)
ws.ReadXML(out=ws.vmr_field_raw, filename=ws.vmrtmp)
# and we clean up the dummy files (not completely, though. but we write an empty
#  variable into them.)
ws.Delete(ws.strtmp)
ws.Touch(ws.strtmp)
ws.WriteXML(ws.tmpformat, ws.strtmp, ws.abstmp, 0)
ws.WriteXML(ws.tmpformat, ws.strtmp, ws.vmrtmp, 0)
