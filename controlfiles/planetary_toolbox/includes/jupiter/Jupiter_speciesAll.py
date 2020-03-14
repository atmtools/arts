#
# ...
#

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
# INCLUDE "parts/Atm-and-SpecAllFromScenario.arts"
# alternatively, add species by hand (e.g., if not ALL in scneario shall be
# considered) and read atm corresponding fields
ws.ArrayOfStringSet(
    ws.aostrtmp,
    ["CO2", "CO", "CH4", "C2H2", "C2H6", "PH3", "H2S", "C2H4", "C3H8", "H2", "He"],
)
ws.execute_controlfile("parts/Atm-and-SpecByHand.arts")
# now appending the species and atm fields, which have other than the basename,
# hence need to be done by hand: H2O, HCN, NH3
### H2O ###
ws.ArrayOfStringSet(ws.aostrtmp, ["H2O"])
ws.StringSet(ws.caseext, ".H2O_high.xml.gz")
ws.execute_controlfile("parts/Atm-and-SpecSingleSpec.arts")
### HCN ###
ws.ArrayOfStringSet(ws.aostrtmp, ["HCN"])
ws.StringSet(ws.caseext, ".HCN_upperlim.xml.gz")
ws.execute_controlfile("parts/Atm-and-SpecSingleSpec.arts")
### NH3 ###
ws.ArrayOfStringSet(ws.aostrtmp, ["NH3"])
ws.StringSet(ws.caseext, ".NH3_high.xml.gz")
ws.execute_controlfile("parts/Atm-and-SpecSingleSpec.arts")
