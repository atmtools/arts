# DEFINITIONS:  -*-sh-*-
#
# Demonstration and test of a transmission calculation in a
# refractive 3D atmosphere including scatterers.
# Airborne type sensor uplooking through a cloud.
#
# Setup is a mixture of TestTransmission and TestRTeCalcMC, adapted for a
# viewing geometry maximizing the differences between the clear and the cloudy
# setup.
#
# 2017-08-30 Jana Mendrok

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
# sensor-only path
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
# Number of Stokes components to be computed
#
ws.IndexSet(ws.stokes_dim, 4)
# Reference ellipsoid
#
ws.refellipsoidEarth(ws.refellipsoid, "WGS84")
# Frequency grid
#
ws.VectorSet(ws.f_grid, array([2.3e11]))
# A pressure grid rougly matching 0 to 80 km, in steps of 2 km.
#
ws.VectorNLogSpace(ws.p_grid, 41, 101300.0, 1.0)
# Atmospheric dimensionality and lat/lon grids
#
# VectorNLinSpace( lat_grid, 11, 5, 13 )
# VectorNLinSpace( lon_grid, 11, -14, -10 )
ws.ReadXML(ws.lat_grid, "artscomponents/montecarlo/lat_grid.xml")
ws.ReadXML(ws.lon_grid, "artscomponents/montecarlo/lon_grid.xml")
ws.AtmosphereSet3D()
# Definition of species
#
ws.abs_speciesSet(species=["H2O-PWR98", "N2-SelfContStandardType", "O2-PWR93"])
# No line data needed here
#
ws.abs_lines_per_speciesSetEmpty()
# Atmospheric profiles
#
ws.AtmRawRead(basename="testdata/tropical")
#
ws.AtmFieldsCalcExpand1D()
# Get ground altitude (z_surface) from z_field
ws.IndexCreate("nlat")
ws.nelemGet(ws.nlat, ws.lat_grid)
ws.IndexCreate("nlon")
ws.nelemGet(ws.nlon, ws.lon_grid)
ws.MatrixSetConstant(ws.z_surface, ws.nlat, ws.nlon, 0.0)
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
# refracted path
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__RefractedPath)
ws.Copy(ws.refr_index_air_agenda, ws.refr_index_air_agenda__GasMicrowavesEarth)
ws.NumericSet(ws.ppath_lmax, 2000.0)
ws.NumericSet(ws.ppath_lraytrace, 500.0)
# Postion and line-of-sight of sensor
#
# VectorSet( rte_pos, [ 0, 5.1, -13.82 ] )
# VectorSet( rte_los, [ 80, 24 ] )
ws.rte_losSet(ws.rte_los, ws.atmosphere_dim, 70.0, 180.0)
ws.rte_posSet(ws.rte_pos, ws.atmosphere_dim, 8000.0, 0.0, 0.0)
ws.VectorSet(ws.rte_pos2, array([], dtype=float64))
# No transmitter position defined
# Radiative transfer agendas
#
ws.Copy(ws.iy_transmitter_agenda, ws.iy_transmitter_agenda__UnitUnpolIntensity)
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Transmission)
# Auxiliary variables
#
# ArrayOfStringSet( iy_aux_vars,
#    [ "Pressure", "Temperature", "VMR, species 0", "Absorption, summed",
#      "Absorption, species 1", "iy", "Transmission", "Optical depth" ] )
ws.abs_xsec_agenda_checkedCalc()
ws.propmat_clearsky_agenda_checkedCalc()
ws.lbl_checkedCalc()
# Set cloud
ws.cloudboxSetManually(
    p1=71617.7922264, p2=17111.6808705, lat1=-1.9, lat2=1.9, lon1=-1.9, lon2=1.9
)
ws.ScatSpeciesInit()
ws.ScatElementsPndAndScatAdd(
    scat_data_files=[
        "testdata/scatData/azi-random_f229-231T214-225r100NP-1ar1_5ice.xml"
    ],
    pnd_field_files=[""],
)
ws.scat_dataCheck()
ws.ScatSpeciesExtendTemperature(T_low=200.0)
ws.scat_dataCalc()
ws.scat_data_checkedCalc()
ws.ReadXML(ws.pnd_field_raw, "artscomponents/montecarlo/pnd_field_raw.xml")
ws.pnd_fieldCalcFrompnd_field_raw()
# Calculate
#
ws.iyCalc()
# WriteXML( output_file_format, iy, "cloudyREFERENCE.xml" )
ws.MatrixCreate("iyREFERENCE")
ws.ReadXML(ws.iyREFERENCE, "cloudyREFERENCE.xml")
ws.Compare(ws.iy, ws.iyREFERENCE, 1e-06)
