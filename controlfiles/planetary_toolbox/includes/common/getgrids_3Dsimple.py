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
# Furthermore, latitude and longitude grids are set to contain the user        #
# specified latitude and longitude limits lat_min, lat_max, lon_min, and       #
# lon_max, respectively. That is, lat_grid and lon_grid only consist of two    #
# points each, sufficient to characterise a homogeneous 3D atmopshere.         #
#                                                                              #
# This file expects the following input parameters:                            #
#   pmin       (Numeric)   Lower limit of pressure to be kept in p_grid.       #
#   pmax       (Numeric)   Upper limit of pressure to be kept in p_grid.       #
#   lat_min    (Numeric)   Lower limit of latitude to be covered by lat_grid.  #
#   lat_max    (Numeric)   Upper limit of latitude to be covered by lat_grid.  #
#   lon_min    (Numeric)   Lower limit of longitude to be covered by lon_grid. #
#   lon_max    (Numeric)   Upper limit of longitude to be covered by lon_grid. #
#   nlat       (Index)     Number of grid points in lat_grid.                  #
#   nlon       (Index)     Number of grid points in lon_grid.                  #
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
ws.VectorNLinSpace(ws.lat_grid, ws.nlat, ws.latmin, ws.latmax)
ws.VectorNLinSpace(ws.lon_grid, ws.nlon, ws.lonmin, ws.lonmax)
