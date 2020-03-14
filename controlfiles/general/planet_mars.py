# DEFINITIONS:  -*-sh-*-
#
############
# Mars specific settings
#
############
#
# Authors: Jana Mendrok
#

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
#
# Isotopologue ratios
#
ws.ReadXML(ws.isotopologue_ratios, "planets/Mars/isotopratio_Mars.xml")
#
# Reference ellipsoid (a spherical ellipsoid must be used for 1D)
#
ws.refellipsoidMars(ws.refellipsoid, "Sphere")
#
# Weight of dry air [g/mol]
# (needed for hydrostatic equilibrium calculations)
# source: http://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
#
ws.NumericSet(ws.molarmass_dry_air, 43.34)
#
# Gravity
# (needed for hydrostatic equilibrium calculations)
#
@arts_agenda
def g0_agenda(ws):
    ws.Ignore(ws.lon)
    ws.Ignore(ws.lat)
    ws.g0Mars()


ws.g0_agenda = g0_agenda

#
# Sidereal rotation period (1.025957 Earth day)
#
ws.NumericSet(ws.planet_rotation_period, 88643.0)
