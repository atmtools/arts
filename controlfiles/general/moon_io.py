# DEFINITIONS:  -*-sh-*-
#
############
# Jupiter specific settings
#
############
#
# Authors: Richard Larsson
#

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
#
# Isotopologue ratios
#
ws.ReadXML(ws.isotopologue_ratios, "planets/Jupiter/isotopratio_Jupiter.xml")
# Note that I do not know this value so I use the Jupiter model as an approximate
#
# Reference ellipsoid (a spherical ellipsoid must be used for 1D)
#
ws.refellipsoidIo(ws.refellipsoid, "Sphere")
#
# Weight of dry air [g/mol]
# (needed for hydrostatic equilibrium calculations)
# source: Number 95% SO2 and 5% SO from Jupiter approximate
#
ws.NumericSet(ws.molarmass_dry_air, 63.110068828)
#
# Gravity
# (needed for hydrostatic equilibrium calculations)
#
@arts_agenda
def g0_agenda(ws):
    ws.Ignore(ws.lon)
    ws.Ignore(ws.lat)
    ws.g0Io()


ws.g0_agenda = g0_agenda

#
# Tidally locked
#
ws.NumericSet(ws.planet_rotation_period, 152853.0)
