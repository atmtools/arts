#
# basic settings for simulations of
#  - passive sensor
#  - 1D
#  - scalar
#  - monochromatic
#  - pencilbeam
#  - clearsky (no clouds/scattering)
#  - radiances
#
#

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# basic settings
################
# scalar RT
ws.IndexSet(ws.stokes_dim, 1)
# 1-dimensional atmosphere
ws.AtmosphereSet1D()
# No jacobian calculation
ws.jacobianOff()
# Clearsky = No scattering
ws.cloudboxOff()
# Monochromatic pencilbeam measurement
ws.sensorOff()
# agenda settings
#################
# (standard) emission calculation
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
# cosmic background radiation
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
# standard surface agenda (i.e., make use of surface_rtprop_agenda)
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
# path defined by sensor line-of-sight
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
