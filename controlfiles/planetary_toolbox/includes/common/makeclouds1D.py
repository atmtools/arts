################################################################################
#                                                                              #
# DO NOT MODIFY this file (unless you are sure what you are doing).            #
# This is only a helper file!                                                  #
#                                                                              #
################################################################################
#                                                                              #
# This file interpolates raw 1D basic atmospheric data (z_field_raw,           #
# t_field_raw, vmr_field_raw) to the calculation grids (p_grid) for 1D         #
# atmosphere output.                                                           #
#                                                                              #
# This file expects the following input parameters:                            #
#   poptmp             (String)        Name of temporary file containing the   #
#                                       list of single scattering data files.  #
#   pndtmp             (String)        Name of temporary file containing       #
#                                       pnd_field_raw.                         #
#   tmpformat          (String)        File format of temporary files.         #
#   atmosphere_dim     as the WSV                                              #
#   f_grid             as the WSV                                              #
#   p_grid             as the WSV                                              #
#                                                                              #
# Output:                                                                      #
#   pnd_field          as the WSV                                              #
#   scat_data          as the WSV                                              #
#   cloudbox_on        as the WSV                                              #
#   cloudbox_limits    as the WSV                                              #
#                                                                              #
################################################################################

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# forth, we get the pnd field and and single scattering data in their respective
#  WSVs pnd_field_raw and scat_data.
ws.ScatSpeciesInit()
ws.ReadXML(out=ws.ssdfiles, filename=ws.poptmp)
ws.ScatSpeciesPndAndScatAdd(scat_data_files=ws.ssdfiles, pnd_fieldarray_file=ws.pndtmp)
ws.scat_dataCheck(scat_data=ws.scat_data_raw)
ws.nelemGet(ws.itmp, ws.p_grid)
ws.IndexStepDown(ws.itmp, ws.itmp)
ws.Extract(ws.pmin_cb, ws.p_grid, ws.itmp)
# prepare a temporary pnd_field over whole atmosphere (very big cloudbox)
ws.cloudboxSetManually(p1=ws.pmax, p2=ws.pmin, lat1=0.0, lat2=0.0, lon1=0.0, lon2=0.0)
ws.pnd_fieldCalcFrompnd_field_raw(zeropadding=1)
# from the whole-atmo pnd_field now determine the actual necessary cloudbox
ws.cloudboxSetAutomatically(particle_field=ws.pnd_field)
ws.pnd_fieldCalcFrompnd_field_raw(zeropadding=1)
# finally, we clean up the dummy files (not completely, though, as deleting
#  files from within ARTS is not possible. instead, we write an empty variable
#  into them.)
ws.Delete(ws.strtmp)
ws.Touch(ws.strtmp)
ws.WriteXML(ws.tmpformat, ws.strtmp, ws.pndtmp, 0)
ws.WriteXML(ws.tmpformat, ws.strtmp, ws.poptmp, 0)
