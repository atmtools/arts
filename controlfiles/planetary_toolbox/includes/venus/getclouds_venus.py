################################################################################
#                                                                              #
# Unless further variables or options for existing variables are introduced,   #
# DO NOT MODIFY this file! This is only a helper file!                         #
#                                                                              #
################################################################################
#                                                                              #
# This file does the actual work of selecting and reading in the RAW cloud     #
# data for Venus as specified by the user. For user specification use, e.g.,   #
# DemoVenusClouds1D.arts as template. The template also contains the detailed  #
# information on which indices are linked to which specific value/selection    #
# for each of the variables. The full arrays, which the indices refer to and   #
# from which the actual values are extracted, are defined in clouds_venus.arts #
# (hence, clouds_venus.arts needs to be included before the file at hand).     #
#                                                                              #
# This file expects the following input parameters:                            #
#   pndcase        (Index)           The cloud field scenario.                 #
#   cloudtypes     (ArrayOfIndex)    The cloud layers to be considered.        #
#                                                                              #
# Files to be included before this file:                                       #
#   includes/venus/clouds_venus.arts                                           #
#   includes/common/createvars.arts                                            #
#                                                                              #
# It provides no direct output, but stores the output in temporary files       #
# (which will the be processed into their final variable containers by         #
# makeclouds1D.arts). This externally stored output is:                        #
#   pnd_field_raw                    as the WSV                                #
#   ssdfiles       (ArrayOfString)   List of files containing the single       #
#                                     scattering data for each particle type.  #
#                                     The list serves as input parameter       #
#                                     scat_data_files to WSM                   #
#                                     ScatSpeciesPndAndScatAdd.                #
# Names of the temporary files are passed into makeclouds1D.arts by the String #
# variables                                                                    #
#   pndtmp                           for pnd_field_raw file                    #
#   poptmp                           for ssdfiles                              #
#                                                                              #
################################################################################

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# We will need to dummy-store some data in files to be able to export data from
# forloops. So we create some dummy names.
# StringSet( tmpformat, "ascii" )
ws.StringSet(ws.tmpformat, "binary")
ws.StringSet(ws.pndtmp, "tmp1.xml")
ws.StringSet(ws.poptmp, "tmp2.xml")
ws.ArrayOfIndexCreate("typesinlayerarray")
# Loop agenda to get pnd & ssd data of all modes in a cloud layer
ws.AgendaCreate("modeloop_agenda")


@arts_agenda
def modeloop_agenda(ws):
    ws.ReadXML(out=ws.pnd_field_raw, filename=ws.pndtmp)
    ws.ReadXML(out=ws.ssdfiles, filename=ws.poptmp)
    ws.Extract(ws.itmp, ws.typesinlayerarray, ws.forloop_index)
    ws.Extract(ws.strtmp, ws.cloudtypearray, ws.itmp)
    ws.Append(ws.pndstr, ws.strtmp)
    ws.Append(ws.ssdstr, ws.strtmp)
    ws.Extract(ws.strtmp, ws.pndarray, ws.pndcase)
    ws.Append(ws.pndstr, ws.strtmp)
    ws.Print(ws.pndstr, 0)
    ws.ReadXML(ws.aogf3tmp, ws.pndstr)
    ws.Append(ws.pnd_field_raw, ws.aogf3tmp)
    ws.WriteXML(ws.tmpformat, ws.pnd_field_raw, ws.pndtmp, 0)
    ws.Append(ws.ssdstr, ws.ssdpostname)
    ws.Print(ws.ssdstr, 0)
    ws.Append(ws.ssdfiles, ws.ssdstr)
    ws.WriteXML(ws.tmpformat, ws.ssdfiles, ws.poptmp, 0)


ws.modeloop_agenda = modeloop_agenda

# Loop agenda to get pnd & ssd data per cloud layer
ws.AgendaCreate("layerloop_agenda")


@arts_agenda
def layerloop_agenda(ws):
    ws.Extract(ws.itmp, ws.cloudtypes, ws.forloop_index)
    ws.Extract(ws.typesinlayerarray, ws.layer2typearray, ws.itmp)
    ws.nelemGet(ws.ncases, ws.typesinlayerarray)
    ws.IndexStepDown(ws.ncases, ws.ncases)
    ws.Copy(ws.forloop_agenda, ws.modeloop_agenda)
    ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)


ws.layerloop_agenda = layerloop_agenda

# Get the cloud pnd field and single scattering data
# ---
# first, create the casename string down to the common filename part in the
# scenario folder. Data is located in:
# Venus.atmo/
ws.Copy(ws.pndstr, ws.cloudbase)
ws.Extract(ws.subatmo, ws.atmoarray, ws.atmo)
ws.Append(ws.pndstr, ws.subatmo)
ws.StringSet(ws.strtmp, "/")
ws.Append(ws.pndstr, ws.strtmp)
ws.Append(ws.pndstr, ws.subatmo)
ws.StringSet(ws.strtmp, ".")
ws.Append(ws.pndstr, ws.strtmp)
ws.Append(ws.pndstr, ws.pndprename)
ws.Copy(ws.ssdstr, ws.cloudbase)
ws.Append(ws.ssdstr, ws.ssdprename)
ws.StringSet(ws.infostr, "Cloud profile data (pnd_field) taken from: ")
ws.Append(ws.infostr, ws.pndstr)
ws.StringSet(ws.strtmp, "*")
ws.Append(ws.infostr, ws.strtmp)
ws.Print(ws.infostr)
# second, we construct the name for the specific data files one-by-one, read the
#  pnd field data and construct pnd_field_raw (array of pnd fields), and store
#  the ssd file names as array in a file. This in order to later use
#  ScatSpeciesPndAndScatAdd to actually get the respective data into the pnd_field_raw
#  and scat_data ARTS WSV. Compared to directly putting the individual pnd
#  field and single scattering data array entries into the WSV one-by-one, this
#  is advantageous as (1) we do not need to read & write the scat_data
#  over and over to get the data out of the ForLoop (saves time) and (2)
#  some basic checks (covering f_grid, correct atmospheric dimension) will be
#  performed on scat_data and pnd_field_raw.
ws.Touch(ws.pnd_field_raw)
ws.Touch(ws.ssdfiles)
ws.WriteXML(ws.tmpformat, ws.pnd_field_raw, ws.pndtmp, 0)
ws.WriteXML(ws.tmpformat, ws.ssdfiles, ws.poptmp, 0)
# third, we loop over the included cloud layers/cloud types.
ws.nelemGet(ws.ncases, ws.cloudtypes)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.Copy(ws.forloop_agenda, ws.layerloop_agenda)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
# all the rest is now done from makeclouds*.arts
