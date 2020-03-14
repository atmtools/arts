################################################################################
#                                                                              #
# Unless further variables or options for existing variables are introduced,   #
# DO NOT MODIFY this file! This is only a helper file!                         #
#                                                                              #
################################################################################
#                                                                              #
# This file does the actual work of selecting and reading in the RAW           #
# wind data for Mars as specified by the user. For user specification use,     #
# e.g., DemoMarsAtmo3D.arts (or its 1D equivalent) as template. The template   #
# also contains the detailed information on which indices are linked to which  #
# specific value/selection for each of the variables. The full arrays, which   #
# the indices refer to and from which the actual values are extracted, are     #
# defined in atmo_mars.arts (hence, atmo_mars.arts needs to be included before #
# the file at hand).                                                           #
#                                                                              #
# This file expects the following input parameters:                            #
#   Ls             (Index)           The season of the atmospheric scenario.   #
#   daytime        (Index)           The daytime of the atmospheric scenario.  #
#   dust           (Index)           The dust loading of the atmospheric       #
#                                     scenario.                                #
#   solar          (Index)           The solar activity level of the           #
#                                     atmospheric scenario.                    #
#   vertwind       (ArrayOfIndex)    Vertical wind setup selected.             #
#   NSwind         (ArrayOfIndex)    N-S wind setup selected.                  #
#   EWwind         (ArrayOfIndex)    E-W wind setup selected.                  #
#                                                                              #
# Files to be included before this file:                                       #
#   includes/mars/atmo_mars.arts                                               #
#   includes/common/createvars.arts                                            #
#                                                                              #
# It provides following output:                                                #
#   wind_u_raw     (GriddedField3)   raw version of wind_u_field               #
#   wind_v_raw     (GriddedField3)   raw version of wind_v_field               #
#   wind_w_raw     (GriddedField3)   raw version of wind_w_field               #
#                                                                              #
################################################################################

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# We will need to dummy-store some data in files to be able to export data from
# forloops. So we create some dummy names.
# StringSet( tmpformat, "ascii" )
ws.StringSet(ws.tmpformat, "binary")
ws.StringSet(ws.Btmp, "tmp1.xml")
# Create a dummy file with empty wind data (because regardless of whether there
#  actually IS data, we are going to read the storage dummy at least once. so we
#  need to create an empty version first.)
ws.Touch(ws.gf3tmp)
# this in case gf3tmp hasn't been used before
ws.Delete(ws.gf3tmp)
# this to throw away possible data in gf3tmp (if it was used before)
ws.Touch(ws.gf3tmp)
# this to initialize it again after deleting
# this to have a properly formated file to read after the forloops
ws.WriteXML(ws.tmpformat, ws.gf3tmp, ws.Btmp, 0)
# Get data for one wind component
ws.AgendaCreate("windloop_agenda")


@arts_agenda
def windloop_agenda(ws):
    ws.Extract(ws.strtmp, ws.casearray, ws.forloop_index)
    ws.Append(ws.specfilename, ws.strtmp)
    #  Print( specfilename, 0 )
    ws.ReadXML(ws.gf3tmp, ws.specfilename)
    ws.WriteXML(ws.tmpformat, ws.gf3tmp, ws.Btmp, 0)


ws.windloop_agenda = windloop_agenda

# Read the wind raw data
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
# now we obtain the wind components one-by-one
# (1) Vertical Wind
ws.Copy(ws.specfilename, ws.atmostr)
ws.Select(ws.casearray, ws.vertwindarray, ws.vertwind)
ws.nelemGet(ws.ncases, ws.casearray)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.Copy(ws.forloop_agenda, ws.windloop_agenda)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
ws.ReadXML(out=ws.wind_w_raw, filename=ws.Btmp)
# (2) N-S Wind
ws.Copy(ws.specfilename, ws.atmostr)
ws.Select(ws.casearray, ws.NSwindarray, ws.NSwind)
ws.nelemGet(ws.ncases, ws.casearray)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.Copy(ws.forloop_agenda, ws.windloop_agenda)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
ws.ReadXML(out=ws.wind_v_raw, filename=ws.Btmp)
# (3) E-W Wind
ws.Copy(ws.specfilename, ws.atmostr)
ws.Select(ws.casearray, ws.EWwindarray, ws.EWwind)
ws.nelemGet(ws.ncases, ws.casearray)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.Copy(ws.forloop_agenda, ws.windloop_agenda)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
ws.ReadXML(out=ws.wind_u_raw, filename=ws.Btmp)
# and we clean up the dummy files (not completely, though. but we write an empty
#  variable into them.)
ws.Delete(ws.strtmp)
ws.Touch(ws.strtmp)
ws.WriteXML(ws.tmpformat, ws.strtmp, ws.Btmp, 0)
