################################################################################
#                                                                              #
# Unless further variables or options for existing variables are introduced,   #
# DO NOT MODIFY this file! This is only a helper file!                         #
#                                                                              #
################################################################################
#                                                                              #
# This file sets up a bunch of arrays that contain names (or pieces of names)  #
# of locations and files of Mars atmospheric data. The entries are supposed to #
# be consistent with the descriptions given to the user in the                 #
# DemoAtmoMars*.arts templates.                                                #
#                                                                              #
################################################################################

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.ArrayOfStringCreate("Lsarray")
ws.ArrayOfStringCreate("daytimearray")
ws.ArrayOfStringCreate("dustarray")
ws.ArrayOfStringCreate("solararray")
ws.ArrayOfStringCreate("basespeciesarray")
ws.ArrayOfStringCreate("basespeciesnamesarray")
ws.ArrayOfStringCreate("CH4array")
ws.ArrayOfStringCreate("H2Oarray")
ws.ArrayOfStringCreate("H2Onamesarray")
ws.ArrayOfStringCreate("Nearray")
ws.ArrayOfStringCreate("vertwindarray")
ws.ArrayOfStringCreate("NSwindarray")
ws.ArrayOfStringCreate("EWwindarray")
ws.StringSet(ws.atmobase, "planets/Mars/MPS/")
ws.ArrayOfStringSet(ws.Lsarray, ["Ls0", "Ls90", "Ls180", "Ls270"])
ws.ArrayOfStringSet(ws.daytimearray, ["day", "night"])
ws.ArrayOfStringSet(ws.dustarray, ["dust-low", "dust-medium", "dust-high"])
ws.ArrayOfStringSet(ws.solararray, ["sol-min", "sol-avg", "sol-max"])
ws.ArrayOfStringSet(
    ws.basespeciesarray,
    [
        "CO",
        "CO2",
        "CO2-CIA-CO2-0",
        "CO2-SelfContPWR93, CO2-ForeignContPWR93",
        "H2",
        "H2O2",
        "H2S",
        "HCl",
        "N2",
        "NO2",
        "O",
        "O2",
        "O3",
        "OCS",
        "OH",
        "SO2",
    ],
)
ws.ArrayOfStringSet(
    ws.basespeciesnamesarray,
    [
        "CO.xml",
        "CO2.xml",
        "CO2.xml",
        "CO2.xml",
        "H2.xml",
        "H2O2.xml",
        "H2S.xml",
        "HCl.xml",
        "N2.xml",
        "NO2.xml",
        "O.xml",
        "O2.xml",
        "O3.xml",
        "OCS.xml",
        "OH.xml",
        "SO2.xml",
    ],
)
ws.ArrayOfStringSet(ws.CH4array, ["CH4.xml", "CH4_high.xml"])
ws.ArrayOfStringSet(ws.H2Oarray, ["H2O-162", "H2O"])
ws.ArrayOfStringSet(ws.H2Onamesarray, ["H2O-162.xml", "H2O.xml"])
ws.ArrayOfStringSet(
    ws.Nearray,
    [
        "SZA.0-30.Ne.xml",
        "SZA.30-50.Ne.xml",
        "SZA.50-70.Ne.xml",
        "SZA.70-90.Ne.xml",
        "SZA.120-180.Ne.xml",
    ],
)
ws.ArrayOfStringSet(ws.vertwindarray, ["wind_w.xml"])
ws.ArrayOfStringSet(ws.NSwindarray, ["wind_v.xml"])
ws.ArrayOfStringSet(ws.EWwindarray, ["wind_u.xml"])
