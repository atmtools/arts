# DEFINITIONS:  -*-sh-*-
#
# Loads an atmosphere and basic absorption data. Absorption calculation is set
# to "on-the-fly".
#
# Planet          : Earth
# Frequency range : 0 - 100 GHz
# Dimensionality  : 3D
# Altitude range  : 0 - 1000 km, 250 m steps up to 80 km and 16 km above.
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
ws.VectorCreate("pionopshere")
ws.VectorNLogSpace(ws.pionopshere, 59, 0.1, 1e-59)
ws.Append(ws.p_grid, ws.pionopshere)
# Tempature and z covering the ionosphere
ws.ReadXML(ws.t_field_raw, "testdata/tropical.expanded.t.xml")
ws.ReadXML(ws.z_field_raw, "testdata/tropical.expanded.z.xml")
ws.GriddedField3Create("raw3")
ws.abs_speciesInit()
# N2
ws.abs_speciesAdd(species=["N2-SelfContStandardType"])
ws.ReadXML(ws.raw3, "testdata/tropical.N2.xml")
ws.Append(ws.vmr_field_raw, ws.raw3)
# O2
ws.abs_speciesAdd(species=["O2-PWR93"])
ws.ReadXML(ws.raw3, "testdata/tropical.O2.xml")
ws.Append(ws.vmr_field_raw, ws.raw3)
# H2O
ws.abs_speciesAdd(species=["H2O-PWR98"])
ws.ReadXML(ws.raw3, "testdata/tropical.H2O.xml")
ws.Append(ws.vmr_field_raw, ws.raw3)
# Free electrons
ws.abs_speciesAdd(species=["free_electrons"])
ws.ReadXML(ws.raw3, "testdata/ne_iri_solmax_spring_12UTC_0latlon.xml")
ws.Append(ws.vmr_field_raw, ws.raw3)
# here, we do LTE only so far, but need to initialize NLTE t-field accordingly
ws.Touch(ws.nlte_field_raw)
# Interpolate to p_grid (VMR is "zero padded")
ws.AtmFieldsCalcExpand1D(vmr_zeropadding=1)
# Surface
ws.Extract(ws.z_surface, ws.z_field, 0)
# Apply HSE
ws.Extract(ws.p_hse, ws.p_grid, 0)
ws.NumericSet(ws.z_hse_accuracy, 0.1)
ws.atmfields_checkedCalc(bad_partition_functions_ok=ws.bad_partition_functions_ok)
ws.z_fieldFromHSE()
# Magnetic field, u-component
ws.ReadXML(ws.raw3, "planets/Earth/IGRF/IGRF11_2010_200km-5deg-5deg.B_u.xml.gz")
ws.GriddedFieldLatLonRegrid(ws.raw3, ws.lat_grid, ws.lon_grid, ws.raw3)
ws.GriddedFieldPRegrid(ws.raw3, ws.p_grid, ws.raw3)
ws.FieldFromGriddedField(ws.mag_u_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.raw3)
# Magnetic field, v-component
ws.ReadXML(ws.raw3, "planets/Earth/IGRF/IGRF11_2010_200km-5deg-5deg.B_v.xml.gz")
ws.GriddedFieldLatLonRegrid(ws.raw3, ws.lat_grid, ws.lon_grid, ws.raw3)
ws.GriddedFieldPRegrid(ws.raw3, ws.p_grid, ws.raw3)
ws.FieldFromGriddedField(ws.mag_v_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.raw3)
# Magnetic field, w-component
ws.ReadXML(ws.raw3, "planets/Earth/IGRF/IGRF11_2010_200km-5deg-5deg.B_w.xml.gz")
ws.GriddedFieldLatLonRegrid(ws.raw3, ws.lat_grid, ws.lon_grid, ws.raw3)
ws.GriddedFieldPRegrid(ws.raw3, ws.p_grid, ws.raw3)
ws.FieldFromGriddedField(ws.mag_w_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.raw3)
# No transitions here
ws.abs_lines_per_speciesSetEmpty()
# Absorption agendas
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)


@arts_agenda
def propmat_clearsky_agenda(ws):
    ws.Ignore(ws.rtp_mag)
    ws.Ignore(ws.rtp_los)
    ws.propmat_clearskyInit()
    ws.propmat_clearskyAddOnTheFly()
    ws.propmat_clearskyAddFaraday()


ws.propmat_clearsky_agenda = propmat_clearsky_agenda
