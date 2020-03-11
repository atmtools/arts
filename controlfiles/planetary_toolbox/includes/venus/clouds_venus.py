################################################################################
#                                                                              #
# Unless further variables or options for existing variables are introduced,   #
# DO NOT MODIFY this file! This is only a helper file!                         #
#                                                                              #
################################################################################
#                                                                              #
# This file sets up a bunch of arrays that contain names (or pieces of names)  #
# of locations and files of Venus cloud data. The entries are supposed to be   #
# consistent with the descriptions given to the user in the                    #
# DemoVenusClouds*.arts templates.                                             #
#                                                                              #
################################################################################

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.StringCreate("pndprename")
ws.StringCreate("ssdprename")
ws.StringCreate("ssdpostname")
ws.ArrayOfStringCreate("cloudtypearray")
ws.ArrayOfArrayOfIndexCreate("layer2typearray")
ws.ArrayOfStringCreate("pndarray")
ws.StringSet(ws.cloudbase, "planets/Venus/SAT/")
ws.StringSet(ws.pndprename, "pnd_field__H2SO4__")
ws.StringSet(ws.ssdprename, "Venus.scat_data__H2SO4__")
ws.StringSet(ws.ssdpostname, "RI-std.xml")
ws.ArrayOfStringSet(
    ws.cloudtypearray,
    [
        "LowerHaze-Mode1-bulk__",
        "LowerCloud-Mode1-bulk__",
        "LowerCloud-Mode2-bulk__",
        "LowerCloud-Mode3-bulk__",
        "MiddleCloud-Mode1-bulk__",
        "MiddleCloud-Mode2-bulk__",
        "MiddleCloud-Mode3-bulk__",
        "UpperCloud-Mode1-bulk__",
        "UpperCloud-Mode2-bulk__",
        "UpperHaze-Mode1-bulk__",
    ],
)
# this sorts the entries of cloudtypearray into groups of modes per cloudlayer
ws.Touch(ws.layer2typearray)
ws.ArrayOfIndexSet(ws.aoitmp, [0])
ws.Append(ws.layer2typearray, ws.aoitmp)
ws.ArrayOfIndexSet(ws.aoitmp, [1, 2, 3])
ws.Append(ws.layer2typearray, ws.aoitmp)
ws.ArrayOfIndexSet(ws.aoitmp, [4, 5, 6])
ws.Append(ws.layer2typearray, ws.aoitmp)
ws.ArrayOfIndexSet(ws.aoitmp, [7, 8])
ws.Append(ws.layer2typearray, ws.aoitmp)
ws.ArrayOfIndexSet(ws.aoitmp, [9])
ws.Append(ws.layer2typearray, ws.aoitmp)
ws.ArrayOfStringSet(
    ws.pndarray, ["KH80-Nstd-box-profile.xml", "KH80-Nalt-box-profile.xml"]
)
