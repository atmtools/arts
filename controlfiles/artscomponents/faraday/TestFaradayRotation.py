# DEFINITIONS:  -*-sh-*-
#
# Demonstration and test of calculation including Faraday rotation
#
# The rotation for a signal transmitted along the zenith direction is
# calculated.
#
# 2013-03-20, Patrick Eriksson

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
# Number of Stokes components to be computed
#
ws.IndexSet(ws.stokes_dim, 4)
# Frequency grid
#
ws.VectorNLogSpace(ws.f_grid, 101, 100000000.0, 5000000000.0)
# VectorSet( f_grid, [0.1e9] )
# Dimensionality of the atmosphere
#
ws.AtmosphereSet1D()
# A pressure grid rougly matching 0 to 1000 km, in steps of 2km.
#
ws.VectorNLogSpace(ws.p_grid, 501, 101300.0, 1e-80)
# Tempature and z covering the ionosphere
#
ws.ReadXML(ws.t_field_raw, "testdata/tropical.expanded.t.xml")
ws.ReadXML(ws.z_field_raw, "testdata/tropical.expanded.z.xml")
#
# Include som gas species
#
# (For demonstration only, no absorption set)
ws.abs_speciesInit()
ws.GriddedField3Create("raw3")
# N2
#
ws.abs_speciesAdd(species=["N2"])
ws.ReadXML(ws.raw3, "testdata/tropical.N2.xml")
ws.Append(ws.vmr_field_raw, ws.raw3)
# O2
#
ws.abs_speciesAdd(species=["O2"])
ws.ReadXML(ws.raw3, "testdata/tropical.O2.xml")
ws.Append(ws.vmr_field_raw, ws.raw3)
# H2O
#
ws.abs_speciesAdd(species=["H2O"])
ws.ReadXML(ws.raw3, "testdata/tropical.H2O.xml")
ws.Append(ws.vmr_field_raw, ws.raw3)
#
# Add free electrons
#
ws.abs_speciesAdd(species=["free_electrons"])
ws.ReadXML(ws.raw3, "testdata/ne_iri_solmax_spring_12UTC_0latlon.xml")
ws.Append(ws.vmr_field_raw, ws.raw3)
ws.Touch(ws.nlte_field_raw)
ws.Touch(ws.nlte_vibrational_energies)
# Interpolate to p_grid (VMR is "zero padded")
#
ws.AtmFieldsCalc(vmr_zeropadding=1)
# Surface altitude
ws.MatrixSetConstant(ws.z_surface, 1, 1, 0.0)
# Apply HSE
#
# (this roughly recreates the original altitudes for the input electron density
#  and magnetic field)
#
ws.VectorSet(ws.lat_true, np.array([0.0]))
ws.VectorSet(ws.lon_true, np.array([0.0]))
ws.NumericSet(ws.p_hse, 101300.0)
ws.NumericSet(ws.z_hse_accuracy, 10.0)
#
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.z_fieldFromHSE()
# Magnetic field components
ws.ReadXML(ws.raw3, "testdata/bu_igrf11_2000_0latlon.xml")
ws.GriddedFieldPRegrid(ws.raw3, ws.p_grid, ws.raw3)
ws.FieldFromGriddedField(ws.mag_u_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.raw3)
ws.ReadXML(ws.raw3, "testdata/bv_igrf11_2000_0latlon.xml")
ws.GriddedFieldPRegrid(ws.raw3, ws.p_grid, ws.raw3)
ws.FieldFromGriddedField(ws.mag_v_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.raw3)
ws.ReadXML(ws.raw3, "testdata/bw_igrf11_2000_0latlon.xml")
ws.GriddedFieldPRegrid(ws.raw3, ws.p_grid, ws.raw3)
ws.FieldFromGriddedField(ws.mag_w_field, ws.p_grid, ws.lat_grid, ws.lon_grid, ws.raw3)
#
# Absorption including Faraday rotation
#
# Agenda for scalar gas absorption calculation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
# No line data needed here
#
ws.abs_lines_per_speciesSetEmpty()
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly_Faraday)
# No jacobian calculations
#
ws.jacobianOff()
# No scattering
#
ws.cloudboxOff()
# Check model atmosphere
#
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
# Propagation path agendas and variables
#
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
ws.NumericSet(ws.ppath_lmax, 10000.0)
# Radiative transfer agendas and variables
#
@arts_agenda
def iy_transmitter_agenda(ws):
    ws.Ignore(ws.rtp_pos)
    ws.Ignore(ws.rtp_los)
    ws.iy_transmitterSinglePol()


ws.iy_transmitter_agenda = iy_transmitter_agenda

ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Transmission)
# Sensor/receiver and transmitter
#
ws.MatrixSet(ws.sensor_pos, np.array([[0.0]]))
ws.MatrixSet(ws.sensor_los, np.array([[0.0]]))
ws.MatrixSet(ws.transmitter_pos, [])
# Dummy value
ws.ArrayOfIndexSet(ws.instrument_pol, [5])
#
ws.sensorOff()
ws.sensor_checkedCalc()
# Auxilary variables
#
ws.ArrayOfStringSet(ws.iy_aux_vars, ["Faraday rotation"])
# Temporarily removed, PE 180313
ws.ArrayOfStringSet(ws.iy_aux_vars, [])
# Calculate
#
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.lbl_checkedCalc()
ws.yCalc()
# Extraxt total Faraday rotation
ws.VectorCreate("farrot_total")
# Temporarily rempved, PE 180313
# Extract( farrot_total, y_aux, 0 )
# WriteXML( "ascii", f_grid,      "f.xml"        )
# WriteXML( "ascii", vmr_field,   "vmr.xml"      )
# WriteXML( "ascii", z_field,     "z_field.xml"  )
# WriteXML( "ascii", mag_w_field, "bw_field.xml" )
# WriteXML( "ascii", farrot_total, "farrot.xml"  )
# WriteXML( "ascii", y, "yREFERENCE.xml" )
# WriteXML( "ascii", farrot_total, "farrot_totalREFERENCE.xml" )
# Expected results
#
ws.VectorCreate("yREFERENCE")
ws.VectorCreate("frtotREFERENCE")
#
ws.ReadXML(ws.yREFERENCE, "yREFERENCE.xml")
ws.ReadXML(ws.frtotREFERENCE, "farrot_totalREFERENCE.xml")
# Check
#
ws.Compare(ws.y, ws.yREFERENCE, 0.0001)
# Temporarily rempved, PE 180313
# Compare( farrot_total, frtotREFERENCE, 1e-2 )
