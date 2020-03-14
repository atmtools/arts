################################################################################
#                                                                              #
# DO NOT MODIFY this file (unless you are sure what you are doing).            #
# This is only a helper file!                                                  #
#                                                                              #
################################################################################
#                                                                              #
# This file provides the grids required for a 3D atmosphere (p/lat/lon_grid).  #
# The pressure grid p_grid is dervied from raw atmosphere's z_field_raw and    #
# reduced to the region between the user specified pressure grid limits p_min  #
# and p_max.                                                                   #
# Furthermore, 5-degree-spaced latitude and longitude grids are set up in the  #
# region between the user specified latitude and longitude limits lat_min,     #
# lat_max, lon_min, and lon_max, respectively.                                 #
#                                                                              #
# This file expects the following input parameters:                            #
#   pmin       (Numeric)   Lower limit of pressure to be kept in p_grid.       #
#   pmax       (Numeric)   Upper limit of pressure to be kept in p_grid.       #
#   lat_min    (Numeric)   Lower limit of latitude to be covered by lat_grid.  #
#   lat_max    (Numeric)   Lower limit of latitude to be covered by lat_grid.  #
#   lon_min    (Numeric)   Lower limit of longitude to be covered by lat_grid. #
#   lon_max    (Numeric)   Lower limit of longitude to be covered by lat_grid. #
#                                                                              #
# Output:                                                                      #
#   p_grid                                                                     #
#   lat_grid                                                                   #
#   lon_grid                                                                   #
################################################################################

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# Derive p_grid from atmosphere data (namely from the z_field data)
ws.p_gridFromZRaw(no_negZ=0)
ws.VectorCrop(ws.p_grid, ws.p_grid, ws.pmin, ws.pmax)
# Define lat and lon grids from given min/max values
# (a) first create a global grids with 5deg spacing (that's what Jupiter's
# magfield, which is the only 3D field in the data for now, has)
ws.VectorLinSpace(ws.lat_grid, -90.0, 90.0, 5.0)
ws.VectorLinSpace(ws.lon_grid, -360.0, 360.0, 5.0)
# (b) crop the grids to the region the user wants
ws.VectorCrop(ws.lat_grid, ws.lat_grid, ws.latmin, ws.latmax)
ws.VectorCrop(ws.lon_grid, ws.lon_grid, ws.lonmin, ws.lonmax)
