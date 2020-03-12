# DEFINITIONS:  -*-sh-*-
#
# ARTS control file for testing 1D propagation path calculations.
#
# It is also demonstration of how to use the ForLopp agenda.
#
# 2012-02-17 Patrick Eriksson

import numpy as np
import arts
from arts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/agendas.arts")
# Agenda for scalar gas absorption calculation
ws.Copy(ws.abs_xsec_agenda, ws.abs_xsec_agenda__noCIA)
# sensor-only path
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
# Number of Stokes components to be computed
#
ws.IndexSet(ws.stokes_dim, 1)
# Reference ellipsoid
#
ws.refellipsoidEarth(ws.refellipsoid, "Sphere")
# A pressure grid rougly matching 0 to 80 km.
#
ws.VectorNLogSpace(ws.p_grid, 41, 100000.0, 1.0)
# Atmospheric dimensionality
#
ws.AtmosphereSet1D()
# Water vapour needed if refraction will be calculated
#
ws.abs_speciesSet(species=["H2O"])
# Read a 1D atmospheric case
#
ws.AtmRawRead(basename="testdata/tropical")
ws.AtmFieldsCalc()
# Ground altitude (z_surface)
#
ws.MatrixSetConstant(ws.z_surface, 1, 1, 500.0)
# No jacobian calculations
#
ws.jacobianOff()
# Initializing cloudbox: No scattering
#
ws.cloudboxOff()
# To instead make tests with cloudbox on, activate those:
# FlagOn( cloudbox_on )
# ArrayOfIndexSet( cloudbox_limits, [ 4, 10 ] )
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
ws.VectorSet(ws.rte_pos, np.array([600000.0]))
# VectorSet( rte_los, [ 112.2550477986 ] )  # Angle that just touches for 600e3
ws.VectorSet(ws.rte_los, np.array([113.0]))
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
            [6.0e05],
            [6.0e05],
            [6.0e05],
            [6.0e05],
            [3.0e03],
            [1.0e04],
            [9.0e03],
            [5.0e02],
            [5.0e02],
        ]
    ),
)
ws.MatrixSet(
    ws.sensor_los,
    np.array(
        [[45.0], [95.0], [113.0], [180.0], [0.0], [90.0], [100.0], [45.0], [100.0]]
    ),
)
ws.IndexCreate("ilast")
ws.nrowsGet(ws.ilast, ws.sensor_pos)
ws.IndexStepDown(ws.ilast, ws.ilast)


@arts_agenda
def forloop_agenda(ws):
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
