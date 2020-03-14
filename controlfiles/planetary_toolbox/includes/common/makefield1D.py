################################################################################
#                                                                              #
# DO NOT MODIFY this file (unless you are sure what you are doing).            #
# This is only a helper file!                                                  #
#                                                                              #
################################################################################
#                                                                              #
# This file interpolates raw 1D GriddedField3 data, e.g., wind component       #
# fields to the calculation grids (p_grid) for 1D atmosphere output.           #
#                                                                              #
# This file expects the following input parameters:                            #
#   p_grid             as the WSV                                              #
#   rawfield           (GriddedField3) a raw atmospheric field (note: the      #
#                                       respective field, e.g. wind_w_raw, has #
#                                       to be copied to rawfield BEFORE        #
#                                       including the file at hand.            #
#   interp_order       (Index)         Grid interpolation order                #
#   vmr_zeropad        (Index)         Flag, whether to fill VMR at            #
#                                       non-covered profile regions with zeros #
# Output:                                                                      #
#   finalfield         (Tensor3)       the atmospheric field interpolated to   #
#                                       p_grid, ready for doing RT calcs.      #
#                                                                              #
################################################################################

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.GriddedFieldPRegrid(ws.rawfield, ws.p_grid, ws.rawfield, 1, ws.auxfield_zeropad)
ws.FieldFromGriddedField(
    ws.finalfield, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.rawfield
)
