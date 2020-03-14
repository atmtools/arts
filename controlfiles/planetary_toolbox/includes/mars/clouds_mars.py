################################################################################
#                                                                              #
# Unless further variables or options for existing variables are introduced,   #
# DO NOT MODIFY this file! This is only a helper file!                         #
#                                                                              #
################################################################################
#                                                                              #
# This file sets up a bunch of arrays that contain names (or pieces of names)  #
# of locations and files of Mars cloud data. The entries are supposed to be    #
# consistent with the descriptions given to the user in the                    #
# DemoMarsClouds*.arts templates.                                              #
#                                                                              #
################################################################################

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.StringCreate("pndprename")
ws.StringCreate("ssdprename")
ws.ArrayOfStringCreate("cloudtypearray")
ws.ArrayOfStringCreate("cloudprofilearray")
ws.ArrayOfStringCreate("RIarray")
ws.ArrayOfIndexCreate("dusttypes")
ws.ArrayOfIndexCreate("co2types")
ws.ArrayOfIndexCreate("h2otypes")
ws.StringSet(ws.cloudbase, "planets/Mars/SAT/")
ws.StringSet(ws.pndprename, "pnd_field__")
ws.StringSet(ws.ssdprename, "Mars.scat_data__")
ws.ArrayOfStringSet(
    ws.cloudtypearray,
    [
        "Dust__small-size-bulk__",
        "Dust__medium-size-bulk__",
        "Dust__large-size-bulk__",
        "Dust__verylarge-size-bulk__",
        "Dust__medium-size-bulk__",
        "Dust__medium-size-bulk__",
        "Dust__medium-size-bulk__",
        "Dust__medium-size-bulk__",
        "Dust__medium-size-bulk__",
        "CO2ice__mesospheric-day-bulk__",
        "CO2ice__polarnight-chan1-bulk__",
        "CO2ice__polarnight-chan1-bulk__",
        "CO2ice__polarnight-chan1-bulk__",
        "CO2ice__polarnight-chan4-bulk__",
        "CO2ice__polarnight-chan4-bulk__",
        "H2Oice__Type1-bulk__",
        "H2Oice__Type1-bulk__",
        "H2Oice__Type2-bulk__",
    ],
)
ws.ArrayOfStringSet(
    ws.cloudprofilearray,
    [
        "verywell-mixed__tau1075cm-1_1e-01.xml",
        "verywell-mixed__tau1075cm-1_1e-01.xml",
        "verywell-mixed__tau1075cm-1_1e-01.xml",
        "verywell-mixed__tau1075cm-1_1e-01.xml",
        "well-mixed__tau1075cm-1_1e-01.xml",
        "mixed__tau1075cm-1_1e-01.xml",
        "confined__tau1075cm-1_1e-01.xml",
        "very-confined__tau1075cm-1_1e-01.xml",
        "highly-confined__tau1075cm-1_1e-01.xml",
        "meso-gauss-profile__tau1024nm_2e-01.xml",
        "20km-narrow-box-profile__ext1024nm_3e-01.xml",
        "8km-narrow-box__ext1024nm_3e-01.xml",
        "8km-tower-box-profile__ext1024nm_3e-01.xml",
        "10km-narrow-box-profile__tau1024nm_1e+00.xml",
        "4km-narrow-box__tau1024nm_1e+00.xml",
        "igh-wide-profile__tau825cm-1_5e-02.xml",
        "low-narrow-profile__tau825cm-1_5e-02.xml",
        "high-wide-profile__tau825cm-1_2e-01.xml",
    ],
)
# this links the cloud scenario arrays to the clopudtypearray above
ws.ArrayOfIndexSet(ws.dusttypes, [0, 1, 2, 3, 4, 5, 6, 7, 8])
ws.ArrayOfIndexSet(ws.co2types, [9, 10, 11, 12, 13, 14])
ws.ArrayOfIndexSet(ws.h2otypes, [15, 16, 17])
ws.ArrayOfStringSet(
    ws.RIarray, ["RI-minabs.xml", "RI-std.xml", "RI-maxabs.xml", "RI-maetzler06.xml"]
)
