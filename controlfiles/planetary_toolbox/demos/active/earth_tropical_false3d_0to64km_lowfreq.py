# DEFINITIONS:  -*-sh-*-
#
# Loads an atmosphere and basic absorption data. Absorption calculation is set
# to "on-the-fly".
#
# Planet          : Earth
# Frequency range : 0 - 100 GHz
# Dimensionality  : 3D
# Altitude range  : 0 - 80 km, 250 m steps
# Gas dataset     : Fascod tropical, expanded to 3D
# Gas species     : H20, O2 and N2, with standard absorption models
#
# Author: Patrick Eriksson

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/planet_earth.arts")
ws.VectorNLogSpace(ws.p_grid, 321, 101300.0, 1.0)
ws.abs_speciesSet(species=["H2O-PWR98", "N2-SelfContStandardType", "O2-PWR93"])
ws.AtmRawRead(basename="testdata/tropical")
ws.AtmFieldsCalcExpand1D()
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.Extract(ws.p_hse, ws.p_grid, 0)
ws.NumericSet(ws.z_hse_accuracy, 0.1)
ws.atmfields_checkedCalc(bad_partition_functions_ok=ws.bad_partition_functions_ok)
ws.z_fieldFromHSE()
ws.abs_lines_per_speciesSetEmpty()
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
