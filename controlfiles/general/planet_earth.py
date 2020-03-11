# DEFINITIONS:  -*-sh-*-
#
############
# Earth specific settings
#
############
#
# Authors: Jana Mendrok
#

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
#
# Isotopologue ratios
#
ws.isotopologue_ratiosInitFromBuiltin()
#
# Reference ellipsoid (a spherical ellipsoid must be used for 1D)
#
ws.refellipsoidEarth(ws.refellipsoid, "Sphere")
#
# Weight of dry air
# (needed for hydrostatic equilibrium calculations)
#
ws.NumericSet(ws.molarmass_dry_air, 28.966)
#
# Gravity
# (needed for hydrostatic equilibrium calculations)
#
@arts_agenda
def g0_agenda(ws):
    ws.Ignore(ws.lon)
    ws.g0Earth()


ws.g0_agenda = g0_agenda

#
# Sidereal rotation period (23h 56min 4.1 s)
#
ws.NumericSet(ws.planet_rotation_period, 86164.1)
