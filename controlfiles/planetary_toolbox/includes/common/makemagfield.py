################################################################################
#                                                                              #
# DO NOT MODIFY this file (unless you are sure what you are doing).            #
# This is only a helper file!                                                  #
#                                                                              #
################################################################################
#                                                                              #
# This file interpolates raw 3D magnetic field data (all three components) to  #
# the calculation grids (p/lat/lon_grid) for 1D or 3D atmosphere output.       #
#                                                                              #
# This file expects the following input parameters:                            #
#   p_grid             as the WSV                                              #
#   lat_true           as the WSV                                              #
#   lon_true           as the WSV                                              #
#   mag_u/v/w_raw      (GriddedField3) raw versions of mag_u/v/w_field         #
#   interp_order       (Index)         Grid interpolation order                #
#   auxfield_zeropad   (Index)         Flag, whether to fill magfield at       #
#                                       non-covered profile regions with zeros #
#                                                                              #
# Output:                                                                      #
#   mag_u/v/w_field    as the WSV                                              #
#                                                                              #
################################################################################

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# do the interpolation/regridding separately for each magfield component
ws.GriddedFieldLatLonRegrid(ws.mag_u_raw, ws.lat_true, ws.lon_true, ws.mag_u_raw, 1)
ws.GriddedFieldPRegrid(ws.mag_u_raw, ws.p_grid, ws.mag_u_raw, 1, ws.auxfield_zeropad)
ws.FieldFromGriddedField(
    ws.mag_u_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.mag_u_raw
)
ws.GriddedFieldLatLonRegrid(ws.mag_v_raw, ws.lat_true, ws.lon_true, ws.mag_v_raw, 1)
ws.GriddedFieldPRegrid(ws.mag_v_raw, ws.p_grid, ws.mag_v_raw, 1, ws.auxfield_zeropad)
ws.FieldFromGriddedField(
    ws.mag_v_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.mag_v_raw
)
ws.GriddedFieldLatLonRegrid(ws.mag_w_raw, ws.lat_true, ws.lon_true, ws.mag_w_raw, 1)
ws.GriddedFieldPRegrid(ws.mag_w_raw, ws.p_grid, ws.mag_w_raw, 1, ws.auxfield_zeropad)
ws.FieldFromGriddedField(
    ws.mag_w_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.mag_w_raw
)
