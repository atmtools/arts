################################################################################
#                                                                              #
# This is a demo/template file. The USER is supposed to MODIFY it according    #
# to his/her needs (better, make a copy of it and adapt the copy).             #
#                                                                              #
################################################################################
#                                                                              #
# This is a template file for deriving Venus cloud data from the arts-xml-data #
# package and convert it to the common spatial grids (p_grid), such that they  #
# can be applied in radiative transfer calculations. It is for a 1D atmosphere.#
#                                                                              #
# It provides following output:                                                #
#   pnd_field         as the WSV                                               #
#   scat_data   as the WSV                                               #
#   cloudbox_on       as the WSV                                               #
#   cloudbox_limits   as the WSV                                               #
#                                                                              #
# The user is supposed to select (cloud layers and density) from lists.        #
# Details of setting rules are given at the place of the settings.             #
#                                                                              #
# Selections and settings to be done are between the flags START USER SETTINGS #
# and END USER SETTINGS. The rest of the file shall not be modified,           #
# particularly settings marked with 'do NOT modify'.                           #
#                                                                              #
# This template makes use of the following include files                       #
#   includes/venus/cloudsatmo_venus.arts                                       #
#   includes/venus/getclouds_venus.arts                                        #
#   includes/common/makeclouds1D.arts                                          #
#                                                                              #
################################################################################

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# do NOT modify
# set up name arrays and the like for selections
ws.execute_controlfile("planetary_toolbox/includes/venus/clouds_venus.arts")
# NOT modify
# prepare the variables for the cloud case (pnd & ssd) selections
ws.IndexCreate("pndcase")
ws.ArrayOfIndexCreate("cloudtypes")
################################################################################
# START USER SETTINGS - Modify selections according to you wishes              #
################################################################################
# Select the basic cloud field scenario to be used
# ---
# For explanation of the cases, see TN3 Sec4.4 and Tab4.2
###
# (0) standard            (KH80, Tab1)
# (1) alternative         (KH80, Tab4 average)
ws.IndexSet(ws.pndcase, 0)
# Select the cloud layers to be included
# ---
# We recommend to use all of them. But you might select a subset, e.g., for
#  sensitivity studies. However, for each of them the full set of particle modes
#  available will be used (if you want to use individual modes only, you can not
#  use this template, but have to make your very own, tailored controlfile (you
#  might ask for help via arts-users mailing list).
###
# lower   lower   middle   upper   upper
#  haze,  cloud,   cloud,  cloud,   haze
#   0  ,    1  ,     2  ,    3  ,    4
# select ALL the ones you want
ws.ArrayOfIndexSet(ws.cloudtypes, [0, 1, 2, 3, 4])
################################################################################
# END USER SETTINGS                                                            #
################################################################################
# do NOT modify
# now, let the prepared include files do the actual work:
# (a) prepare the pnd field and single scattering data of all particle types/
#      cloud layers to include for digestion into respective ARTS WSVs
ws.execute_controlfile("planetary_toolbox/includes/venus/getclouds_venus.arts")
# (b) get the data into the actual WSVs and do conversion of pnd field from raw
#      data with individual grids to the common p/lat/lon_grids
ws.execute_controlfile("planetary_toolbox/includes/common/makeclouds1D.arts")
