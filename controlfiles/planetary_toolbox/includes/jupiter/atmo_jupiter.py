################################################################################
#                                                                              #
# Unless further variables or options for existing variables are introduced,   #
# DO NOT MODIFY this file! This is only a helper file!                         #
#                                                                              #
################################################################################
#                                                                              #
# This file sets up a bunch of arrays that contain names (or pieces of names)  #
# of locations and files of Jupiter atmospheric data. The entries are supposed #
# to be consistent with the descriptions given to the user in the              #
# DemoAtmoJupiter*.arts templates.                                             #
#                                                                              #
################################################################################

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.ArrayOfStringCreate("atmoarray")
ws.ArrayOfStringCreate("basespeciesarray")
ws.ArrayOfStringCreate("basespeciesnamesarray")
ws.ArrayOfStringCreate("H2Oarray")
ws.ArrayOfStringCreate("NH3array")
ws.ArrayOfStringCreate("CH4array")
ws.ArrayOfStringCreate("H2array")
ws.ArrayOfStringCreate("CH4namesarray")
ws.ArrayOfStringCreate("H2namesarray")
ws.ArrayOfStringCreate("Nearray")
ws.ArrayOfStringCreate("windarray")
ws.ArrayOfStringCreate("Barray")
ws.ArrayOfStringCreate("Bcomparray")
ws.StringSet(ws.atmobase, "planets/Jupiter/MPS/")
ws.ArrayOfStringSet(ws.atmoarray, ["Jupiter.mean", "Jupiter.oval"])
ws.ArrayOfStringSet(
    ws.basespeciesarray,
    [
        "C2H2",
        "C2H4",
        "C2H6",
        "C3H8",
        "CO",
        "CO2",
        "H2S",
        "HCN",
        "He",
        "PH3",
        "H2-CIA-H2-0",
        "H2-CIA-He-0",
        "H2-CIA-CH4-0",
        "CH4-CIA-CH4-0",
    ],
)
ws.ArrayOfStringSet(
    ws.basespeciesnamesarray,
    [
        "C2H2.xml.gz",
        "C2H4.xml.gz",
        "C2H6.xml.gz",
        "C3H8.xml.gz",
        "CO.xml.gz",
        "CO2.xml.gz",
        "H2S.xml.gz",
        "HCN_upperlim.xml.gz",
        "He.xml.gz",
        "PH3.xml.gz",
        "H2.xml.gz",
        "H2.xml.gz",
        "H2.xml.gz",
        "CH4.xml.gz",
    ],
)
ws.ArrayOfStringSet(ws.H2Oarray, ["H2O_low.xml.gz", "H2O_high.xml.gz"])
ws.ArrayOfStringSet(ws.NH3array, ["NH3_low.xml.gz", "NH3_high.xml.gz"])
ws.ArrayOfStringSet(ws.CH4array, ["CH4-212", "CH4"])
ws.ArrayOfStringSet(ws.CH4namesarray, ["CH4-212.xml", "CH4.xml.gz"])
ws.ArrayOfStringSet(ws.H2array, ["H2-12", "H2"])
ws.ArrayOfStringSet(ws.H2namesarray, ["H2-12.xml", "H2.xml.gz"])
ws.ArrayOfStringSet(ws.Nearray, ["Ne.low.xml.gz", "Ne.med.xml.gz", "Ne.high.xml.gz"])
ws.ArrayOfStringSet(ws.windarray, ["wind_u.xml.gz", "wind_u.into-thermal.xml.gz"])
ws.ArrayOfStringSet(ws.Barray, ["planets/Jupiter/Khurana/Khurana."])
ws.ArrayOfStringSet(ws.Bcomparray, ["B_w.xml.gz", "B_v.xml.gz", "B_u.xml.gz"])
