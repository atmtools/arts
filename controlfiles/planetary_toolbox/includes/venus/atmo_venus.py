################################################################################
#                                                                              #
# Unless further variables or options for existing variables are introduced,   #
# DO NOT MODIFY this file! This is only a helper file!                         #
#                                                                              #
################################################################################
#                                                                              #
# This file sets up a bunch of arrays that contain names (or pieces of names)  #
# of locations and files of Venus atmospheric data. The entries are supposed   #
# to be consistent with the descriptions given to the user in the              #
# DemoAtmoVenus*.arts templates.                                               #
#                                                                              #
################################################################################

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.ArrayOfStringCreate("atmoarray")
ws.ArrayOfStringCreate("basespeciesarray")
ws.ArrayOfStringCreate("basespeciesnamesarray")
ws.ArrayOfStringCreate("H2Oarray")
ws.ArrayOfStringCreate("HDOarray")
ws.ArrayOfStringCreate("CH4namesarray")
ws.ArrayOfStringCreate("Nearray")
ws.ArrayOfStringCreate("NSwindarray")
ws.ArrayOfStringCreate("EWwindarray")
ws.StringSet(ws.atmobase, "planets/Venus/MPS/")
ws.ArrayOfStringSet(
    ws.atmoarray,
    [
        "Venus.spicav.night",
        "Venus.spicav.night_cold",
        "Venus.vira.night",
        "Venus.vira.day",
        "Venus.vira.day_highlat",
    ],
)
ws.ArrayOfStringSet(
    ws.basespeciesarray,
    [
        "CO",
        "CO2",
        "CO2-CIA-CO2-0",
        "CO2-SelfContPWR93, CO2-ForeignContPWR93",
        "H2SO4",
        "HCl",
        "HF",
        "N2",
        "NO",
        "NO2",
        "O",
        "O2",
        "O3",
        "OCS",
        "SO",
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
        "H2SO4.xml",
        "HCl.xml",
        "HF.xml",
        "N2.xml",
        "NO.xml",
        "NO2.xml",
        "O.xml",
        "O2.xml",
        "O3.xml",
        "OCS.xml",
        "SO.xml",
        "SO2.xml",
    ],
)
ws.ArrayOfStringSet(ws.H2Oarray, ["H2O_low.xml", "H2O_mid.xml", "H2O_high.xml"])
ws.ArrayOfStringSet(
    ws.HDOarray,
    [
        "H2O-162_low.xml",
        "H2O-162_mid.xml",
        "H2O-162_high.xml",
        "H2O-162_uncorrected.xml",
    ],
)
ws.ArrayOfStringSet(
    ws.Nearray,
    [
        "SZA.0-30.Ne.xml",
        "SZA.30-50.Ne.xml",
        "SZA.50-70.Ne.xml",
        "SZA.70-80.Ne.xml",
        "SZA.80-90.Ne.xml",
        "SZA.90-100.Ne.xml",
        "SZA.100-120.Ne.xml",
    ],
)
ws.ArrayOfStringSet(ws.NSwindarray, ["wind_v.xml"])
ws.ArrayOfStringSet(
    ws.EWwindarray, ["wind_u_min.xml", "wind_u_mid.xml", "wind_u_max.xml"]
)
