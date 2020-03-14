#####
#
# The only task of this file is to create variables that are used by:
#   PLANET/getatmo_PLANET.arts
#   PLANET/getwind_PLANET.arts
#   PLANET/getmagfield_PLANET.arts
#   PLANET/getclouds_PLANET.arts
#   PLANET/atmo_PLANET.arts
#   PLANET/clouds_PLANET.arts
#   common/makefield1D.arts
#   common/makefield3D.arts
#   common/makeclouds1D.arts
#
# It needs to be included ONCE (and ONLY once!), before the above named
# include files are used.
#
#####

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.StringCreate("infostr")
ws.StringCreate("tmpformat")
ws.StringCreate("atmobase")
ws.StringCreate("superatmo")
ws.StringCreate("subatmo")
ws.StringCreate("cloudbase")
ws.StringCreate("atmostr")
ws.StringCreate("specfilename")
ws.ArrayOfStringCreate("speciesname")
ws.ArrayOfIndexCreate("case")
ws.ArrayOfStringCreate("casearray")
ws.IndexCreate("ncases")
ws.NumericCreate("pmin")
ws.NumericCreate("pmax")
ws.NumericCreate("pmax_cb")
ws.NumericCreate("pmin_cb")
ws.StringCreate("pndstr")
ws.StringCreate("ssdstr")
ws.ArrayOfStringCreate("ssdfiles")
ws.StringCreate("strtmp")
ws.StringCreate("Btmp")
ws.StringCreate("abstmp")
ws.StringCreate("vmrtmp")
ws.StringCreate("pndtmp")
ws.StringCreate("poptmp")
ws.IndexCreate("itmp")
ws.ArrayOfIndexCreate("aoitmp")
ws.GriddedField3Create("gf3tmp")
ws.ArrayOfGriddedField3Create("aogf3tmp")
ws.IndexCreate("bad_partition_functions_ok")
ws.GriddedField3Create("wind_u_raw")
ws.GriddedField3Create("wind_v_raw")
ws.GriddedField3Create("wind_w_raw")
ws.GriddedField3Create("mag_u_raw")
ws.GriddedField3Create("mag_v_raw")
ws.GriddedField3Create("mag_w_raw")
ws.GriddedField3Create("rawfield")
ws.Tensor3Create("finalfield")
