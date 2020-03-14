################################################################################
#                                                                              #
# Unless further variables or options for existing variables are introduced,   #
# DO NOT MODIFY this file! This is only a helper file!                         #
#                                                                              #
################################################################################
#                                                                              #
# This file does the actual work of selecting and reading in the RAW           #
# magnetic field data for Jupiter as specified by the user. For user           #
# specification use, e.g., DemoJupiterAtmo3D.arts (or its 1D equivalent) as    #
# template. The template also contains the detailed information on which       #
# indices are linked to which specific value/selection for each of the         #
# variables. The full arrays, which the indices refer to and from which the    #
# actual values are extracted, are defined in atmo_mars.arts (hence,           #
# atmo_mars.arts needs to be included before the file at hand).                #
#                                                                              #
# This file expects the following input parameters:                            #
#   atmo           (Index)           The atmospheric scenario.                 #
#   Bcase          (ArrayOfIndex)    Magnetic field setup selected             #
#                                     (off/Khurana).                           #
#                                                                              #
# Files to be included before this file:                                       #
#   includes/jupiter/atmo_jupiter.arts                                         #
#   includes/common/createvars.arts                                            #
#                                                                              #
# It provides following output:                                                #
#   mag_u/v/w_raw  (GriddedField3)   raw versions of mag_u/v/w_field           #
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
# Create a dummy file with empty magfield data (because regardless of whether
#  there actually IS data, we are going to read the storage dummy at least once.
#  so we need to create an empty version first.)
ws.Touch(ws.aogf3tmp)
# this in case aogf3tmp hasn't been used before
ws.Delete(ws.aogf3tmp)
# this to throw away possible data in aogf3tmp (if it was used before)
ws.Touch(ws.gf3tmp)
# this in case gf3tmp hasn't been used before
ws.Delete(ws.gf3tmp)
# this to throw away possible data in gf3tmp (if it was used before)
ws.Touch(ws.gf3tmp)
# this to initialize it again after deleting
ws.Append(ws.aogf3tmp, ws.gf3tmp)
ws.Append(ws.aogf3tmp, ws.gf3tmp)
ws.Append(ws.aogf3tmp, ws.gf3tmp)
# this to have a properly formated file to read after the forloops
ws.WriteXML(ws.tmpformat, ws.aogf3tmp, ws.Btmp, 0)
# Get data for one magfield component
ws.AgendaCreate("Bcomploop_agenda")


@arts_agenda
def Bcomploop_agenda(ws):
    ws.ReadXML(out=ws.aogf3tmp, filename=ws.Btmp)
    ws.Extract(ws.strtmp, ws.Bcomparray, ws.forloop_index)
    ws.Append(ws.specfilename, ws.strtmp)
    #  Print( specfilename, 0 )
    ws.ReadXML(ws.gf3tmp, ws.specfilename)
    ws.Append(ws.aogf3tmp, ws.gf3tmp)
    ws.WriteXML(ws.tmpformat, ws.aogf3tmp, ws.Btmp, 0)


ws.Bcomploop_agenda = Bcomploop_agenda

# Get data for the magfield case, looping over the components
ws.AgendaCreate("Bloop_agenda")


@arts_agenda
def Bloop_agenda(ws):
    ws.Touch(ws.aogf3tmp)
    ws.WriteXML(ws.tmpformat, ws.aogf3tmp, ws.Btmp, 0)
    ws.Extract(ws.specfilename, ws.casearray, ws.forloop_index)
    ws.nelemGet(ws.ncases, ws.Bcomparray)
    ws.IndexStepDown(ws.ncases, ws.ncases)
    ws.Copy(ws.forloop_agenda, ws.Bcomploop_agenda)
    ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)


ws.Bloop_agenda = Bloop_agenda

# Read the atmospheric setup
# ---
# (1) Magnetic Field
ws.Select(ws.casearray, ws.Barray, ws.Bcase)
ws.nelemGet(ws.ncases, ws.casearray)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.Copy(ws.forloop_agenda, ws.Bloop_agenda)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
ws.ReadXML(out=ws.aogf3tmp, filename=ws.Btmp)
ws.Extract(ws.mag_w_raw, ws.aogf3tmp, 3)
ws.Extract(ws.mag_v_raw, ws.aogf3tmp, 4)
ws.Extract(ws.mag_u_raw, ws.aogf3tmp, 5)
# and we clean up the dummy files (not completely, though. but we write an empty
#  variable into them.)
ws.Delete(ws.strtmp)
ws.Touch(ws.strtmp)
ws.WriteXML(ws.tmpformat, ws.strtmp, ws.Btmp, 0)
