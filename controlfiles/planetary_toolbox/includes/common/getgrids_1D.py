################################################################################
#                                                                              #
# DO NOT MODIFY this file (unless you are sure what you are doing).            #
# This is only a helper file!                                                  #
#                                                                              #
################################################################################
#                                                                              #
# This file provides the grids required for a 1D atmosphere (p_grid).          #
# The pressure grid p_grid is dervied from raw atmosphere's z_field_raw and    #
# reduced to the region between the user specified pressure grid limits p_min  #
# and p_max.                                                                   #
# Latitude and longitude grids are set empty as required for a 1D case.        #
#                                                                              #
# This file expects the following input parameters:                            #
#   pmin       (Numeric)   Lower limit of pressure to be kept in p_grid.       #
#   pmax       (Numeric)   Upper limit of pressure to be kept in p_grid.       #
#                                                                              #
# Output:                                                                      #
#   p_grid                                                                     #
#   lat_grid  (set as empty)                                                   #
#   lon_grid  (set as empty)                                                   #
#                                                                              #
################################################################################

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# Derive p_grid from atmosphere data (namely from the z_field data)
ws.p_gridFromZRaw(no_negZ=0)
ws.VectorCrop(ws.p_grid, ws.p_grid, ws.pmin, ws.pmax)
ws.VectorSet(ws.lat_grid, [])
ws.VectorSet(ws.lon_grid, [])
