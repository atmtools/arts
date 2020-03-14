# DEFINITIONS:  -*-sh-*-
#
############
# Jupiter specific settings
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
ws.ReadXML(ws.isotopologue_ratios, "planets/Jupiter/isotopratio_Jupiter.xml")
#
# Reference ellipsoid (a spherical ellipsoid must be used for 1D)
#
ws.refellipsoidJupiter(ws.refellipsoid, "Sphere")
#
# Weight of dry air [g/mol]
# (needed for hydrostatic equilibrium calculations)
# source: http://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html
#
ws.NumericSet(ws.molarmass_dry_air, 2.22)
#
# Gravity
# (needed for hydrostatic equilibrium calculations)
#
@arts_agenda
def g0_agenda(ws):
    ws.Ignore(ws.lon)
    ws.Ignore(ws.lat)
    ws.g0Jupiter()


ws.g0_agenda = g0_agenda

#
# Sidereal rotation period (9 h 55 m 30 s)
#
ws.NumericSet(ws.planet_rotation_period, 35730.0)
