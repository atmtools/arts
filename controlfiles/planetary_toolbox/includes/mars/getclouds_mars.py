################################################################################
#                                                                              #
# Unless further variables or options for existing variables are introduced,   #
# DO NOT MODIFY this file! This is only a helper file!                         #
#                                                                              #
################################################################################
#                                                                              #
# This file does the actual work of selecting and reading in the RAW cloud     #
# data for Mars as specified by the user. For user specification use, e.g.,    #
# DemoMarsClouds1D.arts as template. The template also contains the detailed   #
# information on which indices are linked to which specific value/selection    #
# for each of the variables. The full arrays, which the indices refer to and   #
# from which the actual values are extracted, are defined in clouds_mars.arts  #
# (hence, clouds_mars.arts needs to be included before the file at hand).      #
#                                                                              #
# This file expects the following input parameters:                            #
#   dustcase       (ArrayOfIndex)    Selected dust scenario(s).                #
#   co2case        (ArrayOfIndex)    Selected CO2 ice cloud scenario(s).       #
#   h2ocase        (ArrayOfIndex)    Selected H2O ice cloud scenario(s).       #
#   dustRI         (Index)           Refractive index setup for Martian dust.  #
#   co2RI          (Index)           Refractive index setup for CO2 ice.       #
#   h2oRI          (Index)           Refractive index setup for H2O ice.       #
#                                                                              #
# Files to be included before this file:                                       #
#   includes/mars/clouds_mars.arts                                             #
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
ws.IndexCreate("RIcase")
ws.ArrayOfIndexCreate("cloudcase")
ws.ArrayOfIndexCreate("cloudcasearray")
# Loop agenda to get pnd & ssd data of all modes in a cloud layer
ws.AgendaCreate("caseloop_agenda")


@arts_agenda
def caseloop_agenda(ws):
    ws.ReadXML(out=ws.pnd_field_raw, filename=ws.pndtmp)
    ws.ReadXML(out=ws.ssdfiles, filename=ws.poptmp)
    ws.Extract(ws.itmp, ws.cloudcase, ws.forloop_index)
    ws.Extract(ws.itmp, ws.cloudcasearray, ws.itmp)
    ws.Extract(ws.strtmp, ws.cloudtypearray, ws.itmp)
    ws.Append(ws.pndstr, ws.strtmp)
    ws.Append(ws.ssdstr, ws.strtmp)
    ws.Extract(ws.strtmp, ws.cloudprofilearray, ws.itmp)
    ws.Append(ws.pndstr, ws.strtmp)
    ws.Extract(ws.strtmp, ws.RIarray, ws.RIcase)
    ws.Append(ws.ssdstr, ws.strtmp)
    ws.Print(ws.pndstr, 0)
    ws.ReadXML(ws.aogf3tmp, ws.pndstr)
    ws.Append(ws.pnd_field_raw, ws.aogf3tmp)
    ws.WriteXML(ws.tmpformat, ws.pnd_field_raw, ws.pndtmp, 0)
    ws.Print(ws.ssdstr, 0)
    ws.Append(ws.ssdfiles, ws.ssdstr)
    ws.WriteXML(ws.tmpformat, ws.ssdfiles, ws.poptmp, 0)


ws.caseloop_agenda = caseloop_agenda

# Get the cloud pnd field and single scattering data
# ---
# first, create the casename string down to the common filename part in the
# scenario folder. Data is located in:
# Mars.atmo/
# base name for pndfields
ws.Copy(ws.pndstr, ws.cloudbase)
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
ws.Append(ws.pndstr, ws.atmostr)
ws.StringSet(ws.strtmp, "/")
ws.Append(ws.pndstr, ws.strtmp)
# append base file name (Mars.Ls.daytime.dust.) to path construction
ws.Append(ws.pndstr, ws.atmostr)
ws.StringSet(ws.strtmp, ".")
ws.Append(ws.pndstr, ws.strtmp)
# append pndfield marker
ws.Append(ws.pndstr, ws.pndprename)
# base name for scatdata
ws.Copy(ws.ssdstr, ws.cloudbase)
ws.Append(ws.ssdstr, ws.ssdprename)
ws.StringSet(ws.infostr, "Cloud profile data (pnd_field) taken from: ")
ws.Append(ws.infostr, ws.pndstr)
ws.StringSet(ws.strtmp, "*")
ws.Append(ws.infostr, ws.strtmp)
ws.Print(ws.infostr)
# now, we construct the name for the specific data files one-by-one, read the
#  pnd field data and construct pnd_field_raw (array of pnd fields), and store
#  the ssd file names as array in a file. This in order to later use
#  ScatSpeciesPndAndScatAdd to actually get the respective data into the pnd_field_raw
#  and scat_data ARTS WSV. Compared to directly putting the individual pnd
#  field and single scattering data array entries into the WSV one-by-one, this
#  is advantageous as (1) we do not need to read & write the scat_data
#  over and over to get the data out of the ForLoop (saves time) and (2)
#  some basic checks (covering f_grid, correct atmospheric dimension) will be
#  performed on scat_data and pnd_field_raw.
# second, we need to initialize the containers for the pndfield and scatdata
ws.Touch(ws.pnd_field_raw)
ws.Touch(ws.ssdfiles)
ws.WriteXML(ws.tmpformat, ws.pnd_field_raw, ws.pndtmp, 0)
ws.WriteXML(ws.tmpformat, ws.ssdfiles, ws.poptmp, 0)
# third, we loop over the included cloud layers/cloud types. separately per
#  (super)cloud type. no check done, whether indeed only one scenario per
#  (super)cloud type is defined. it's running, if more are defined, but the
#  constructed case doesn't make sense then.
# first, the dust
ws.Copy(ws.cloudcase, ws.dustcase)
ws.Copy(ws.cloudcasearray, ws.dusttypes)
ws.Copy(ws.RIcase, ws.dustRI)
ws.nelemGet(ws.ncases, ws.cloudcase)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.Copy(ws.forloop_agenda, ws.caseloop_agenda)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
# then, the co2 ice
ws.Copy(ws.cloudcase, ws.co2case)
ws.Copy(ws.cloudcasearray, ws.co2types)
ws.Copy(ws.RIcase, ws.co2RI)
ws.nelemGet(ws.ncases, ws.cloudcase)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.Copy(ws.forloop_agenda, ws.caseloop_agenda)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
# eventually, the h2o ice
ws.Copy(ws.cloudcase, ws.h2ocase)
ws.Copy(ws.cloudcasearray, ws.h2otypes)
ws.Copy(ws.RIcase, ws.h2oRI)
ws.nelemGet(ws.ncases, ws.cloudcase)
ws.IndexStepDown(ws.ncases, ws.ncases)
ws.Copy(ws.forloop_agenda, ws.caseloop_agenda)
ws.ForLoop(ws.forloop_agenda, 0, ws.ncases, 1)
# all the rest is now done by makeclouds*.arts
