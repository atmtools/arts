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
#   p_grid             as the WSV                                              #
#   z_field_raw        as the WSV                                              #
#   t_field_raw        as the WSV                                              #
#   vmr_field_raw      as the WSV                                              #
#   interp_order       (Index)         Grid interpolation order                #
#   vmr_zeropad        (Index)         Flag, whether to fill VMR at            #
#                                       non-covered profile regions with zeros #
# Output:                                                                      #
#   z_field            as the WSV                                              #
#   t_field            as the WSV                                              #
#   vmr_field          as the WSV                                              #
#                                                                              #
################################################################################

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# here, we do LTE only so far, but need to initialize NLTE t-field accordingly
ws.Touch(ws.nlte_field_raw)
# this brings z_, t_, and vmr_field_raw to p_grid
ws.AtmFieldsCalc(interp_order=ws.interp_order, vmr_zeropadding=ws.vmr_zeropad)
