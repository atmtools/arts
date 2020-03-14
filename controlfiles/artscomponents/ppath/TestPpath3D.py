# DEFINITIONS:  -*-sh-*-
#
# ARTS control file for testing 3D propagation path calculations
#
# It is also demonstration of how to use the ForLopp agenda.
#
# 2012-02-17 Patrick Eriksson

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
# Agenda for scalar gas absorption calculation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
# sensor-only path
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
# Number of Stokes components to be computed
#
ws.IndexSet(ws.stokes_dim, 1)
# A pressure grid rougly matching 0 to 80 km.
#
ws.VectorNLogSpace(ws.p_grid, 41, 100000.0, 1.0)
# Atmospheric dimensionality and lat/lon grids
#
ws.VectorNLinSpace(ws.lat_grid, 21, 35.0, 55.0)
ws.VectorNLinSpace(ws.lon_grid, 21, -40.0, 40.0)
ws.AtmosphereSet3D()
# Water vapour needed if refraction will be calculated
#
ws.abs_speciesSet(species=["H2O"])
# Read a 1D atmospheric case and expand to *atmosphere_dim*
#
ws.AtmRawRead(basename="testdata/tropical")
ws.AtmFieldsCalcExpand1D()
# Ground altitude (z_surface)
#
ws.IndexCreate("nlat")
ws.IndexCreate("nlon")
ws.nelemGet(ws.nlat, ws.lat_grid)
ws.nelemGet(ws.nlon, ws.lon_grid)
#
ws.MatrixSetConstant(ws.z_surface, ws.nlat, ws.nlon, 500.0)
# No jacobian calculations
#
ws.jacobianOff()
# Initializing cloudbox: No scattering
#
ws.cloudboxOff()
# Activate to make tests with cloudbox on:
# FlagOn( cloudbox_on )
# ArrayOfIndexSet( cloudbox_limits, [ 4, 7, 10, 11, 8, 12 ] )
# A dummy frequency grid
#
ws.VectorSet(ws.f_grid, np.array([1.0e10]))
# Check if atmosphere OK
#
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
# Max step length for the representation of the propagation path
#
ws.NumericSet(ws.ppath_lmax, 20000.0)
#
# A single propagation path:
#
# Set a observation position and line-of-sight (LOS)
#
ws.VectorSet(ws.rte_pos, np.array([6.00e05, 7.01e01, 0.00e00]))
ws.VectorSet(ws.rte_los, np.array([113.0, 180.0]))
ws.VectorSet(ws.rte_pos2, [])
# No transmitter involved
#
# no refraction
#
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
# Calculate propagation path
#
ws.ppathCalc()
# Print complete ppath
#
# Print( ppath, 0 )
# Uncomment to skip batch part below
#
# Exit()
#
# Run through a number of cases, that should run without error
#
# Use sensor_pos/los to store each case
#
ws.MatrixSet(
    ws.sensor_pos,
    np.array(
        [
            [6.00e05, 2.02e01, 2.30e01],
            [6.00e05, 7.03e01, -1.20e01],
            [6.00e05, 2.01e01, -3.00e01],
            [6.00e05, 7.02e01, 1.10e01],
            [6.00e05, 4.10e01, 1.50e00],
            [6.00e05, 3.74e01, -2.00e01],
            [6.00e05, 5.00e01, 1.00e01],
            [6.00e03, 4.01e01, 3.46e01],
            [6.00e03, 4.01e01, -3.90e01],
            [6.00e03, 4.10e01, 1.50e00],
            [6.00e03, 3.74e01, -2.00e01],
            [6.00e03, 5.00e01, 1.00e01],
            [5.00e02, 4.10e01, 2.30e-01],
            [5.00e02, 3.76e01, -1.10e00],
        ]
    ),
)
ws.MatrixSet(
    ws.sensor_los,
    np.array(
        [
            [45.0, 15.0],
            [95.0, 135.0],
            [113.0, 12.0],
            [113.0, -168.0],
            [180.0, 0.0],
            [180.0, 156.0],
            [180.0, -29.0],
            [34.0, -160.0],
            [156.0, 90.0],
            [0.0, 0.0],
            [0.0, 156.0],
            [0.0, -29.0],
            [10.0, 0.0],
            [112.0, 55.0],
        ]
    ),
)
ws.IndexCreate("ilast")
ws.nrowsGet(ws.ilast, ws.sensor_pos)
ws.IndexStepDown(ws.ilast, ws.ilast)


@arts_agenda
def forloop_agenda(ws):
    # Print( forloop_index, 0 )
    ws.VectorExtractFromMatrix(ws.rte_pos, ws.sensor_pos, ws.forloop_index, "row")
    ws.VectorExtractFromMatrix(ws.rte_los, ws.sensor_los, ws.forloop_index, "row")
    ws.ppathCalc()


ws.forloop_agenda = forloop_agenda

#
# no refraction
#
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
ws.ForLoop(ws.forloop_agenda, 0, ws.ilast, 1)
#
# Repeat with refraction
#
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__RefractedPath)
ws.NumericSet(ws.ppath_lraytrace, 1000.0)
ws.Copy(ws.refr_index_air_agenda, ws.refr_index_air_agenda__GasMicrowavesEarth)
ws.ForLoop(ws.forloop_agenda, 0, ws.ilast, 1)
#
# Repeat all with non-spherical ellipsoid
#
ws.refellipsoidEarth(ws.refellipsoid, "WGS84")
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
ws.ForLoop(ws.forloop_agenda, 0, ws.ilast, 1)
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__RefractedPath)
ws.Copy(ws.refr_index_air_agenda, ws.refr_index_air_agenda__GasMicrowavesEarth)
ws.ForLoop(ws.forloop_agenda, 0, ws.ilast, 1)
