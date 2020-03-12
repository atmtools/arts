# DEFINITIONS:  -*-sh-*-
#
# Demonstration and test of a single radio link calculations, using iyCalc.
#
# The link between a satellite and a surface point is simulated. The refracted
# path is at about 5 degress above the horizon (at the surface point).
#
# 2012-08-21, Patrick Eriksson

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/continua.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
# Agenda for scalar gas absorption calculation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
# on-the-fly absorption
ws.Copy(ws.propmat_clearsky_agenda, ws.propmat_clearsky_agenda__OnTheFly)
# Number of Stokes components to be computed
#
ws.IndexSet(ws.stokes_dim, 1)
# Frequency grid
#
ws.VectorSet(ws.f_grid, np.array([1.0e09, 5.0e09, 1.0e10]))
# A pressure grid rougly matching 0 to 80 km, in steps of 250.
#
ws.VectorNLogSpace(ws.p_grid, 321, 101300.0, 1.0)
# Definition of species
#
ws.abs_speciesSet(species=["H2O-PWR98", "N2-SelfContStandardType", "O2-PWR93"])
# No line data needed here
#
ws.abs_lines_per_speciesSetEmpty()
# Dimensionality of the atmosphere
#
ws.AtmosphereSet1D()
# Atmospheric profiles
#
ws.AtmRawRead(basename="testdata/tropical")
#
ws.AtmFieldsCalc(
    ws.t_field,
    ws.z_field,
    ws.vmr_field,
    ws.nlte_field,
    ws.p_grid,
    ws.lat_grid,
    ws.lon_grid,
    ws.t_field_raw,
    ws.z_field_raw,
    ws.vmr_field_raw,
    ws.t_field_raw,
    ws.atmosphere_dim,
    np.array([3.0]),
)
# Surface altitude
ws.MatrixSetConstant(ws.z_surface, 1, 1, 0.0)
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
# transmitter-receiver path
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__TransmitterReceiverPath)
# refracted path
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__RefractedPath)
# The reference values are calculated with Thayer's values for k1, k2 and k3:
@arts_agenda
def refr_index_air_agenda(ws):
    ws.Ignore(ws.f_grid)
    ws.NumericSet(ws.refr_index_air, 1.0)
    ws.NumericSet(ws.refr_index_air_group, 1.0)
    ws.refr_index_airMicrowavesEarth(k1=7.76e-07, k2=6.48e-07, k3=0.003776)


ws.refr_index_air_agenda = refr_index_air_agenda

#
ws.NumericSet(ws.ppath_lmax, 10000.0)
ws.NumericSet(ws.ppath_lraytrace, 100.0)
# Radiative transfer agendas and variables
#
ws.VectorSet(ws.rte_los, [])
# Dummy value
#
ws.Copy(ws.iy_transmitter_agenda, ws.iy_transmitter_agenda__UnitUnpolIntensity)
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Radiolink)
# Postion of sensor/receiver and transmitter
#
ws.VectorSet(ws.rte_pos, np.array([80000.0]))
ws.VectorSet(ws.rte_pos2, np.array([0.0, 5.1]))
# Auxilary variables
#
ws.ArrayOfStringSet(
    ws.iy_aux_vars,
    [
        "Pressure",
        "Temperature",
        "VMR, species 0",
        "Atmospheric loss",
        "Absorption, summed",
        "Free space loss",
        "Free space attenuation",
        "Defocusing loss",
        "Extra path delay",
        "Bending angle",
    ],
)
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
# Run and save
#
ws.iyCalc()
# Set "Atmospheric loss" is main variable:
# (If yCalc is used, iyReplaceFromAux shall be included in iy_main_agenda)
#
ws.iyReplaceFromAux(
    ws.iy, ws.iy_aux, ws.iy_aux_vars, ws.jacobian_do, "Atmospheric loss"
)
# WriteXML( "ascii", iy, "iyREFERENCE.xml" )
# WriteXML( "ascii", ppath, "ppath.xml" )
# OK?
# ---
ws.MatrixCreate("iyREFERENCE")
ws.ReadXML(ws.iyREFERENCE, "iyREFERENCE.xml")
ws.Compare(ws.iy, ws.iyREFERENCE, 1e-05)
