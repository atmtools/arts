# DEFINITIONS:  -*-sh-*-
#
# Demonstrates amd tests rotation of the Stokes vector as a function of scan
# angle. The cfile also shows how an atmosphere completly lacking attenuation
# can be generated.
#
# Authors: Patrick Eriksson

import numpy as np
import pyarts
from pyarts.workspace import Workspace, arts_agenda

ws = Workspace(verbosity=0)
ws.execute_controlfile("general/general.arts")
ws.execute_controlfile("general/agendas.arts")
ws.execute_controlfile("general/planet_earth.arts")
# cosmic background radiation
ws.Copy(ws.iy_space_agenda, ws.iy_space_agenda__CosmicBackground)
# standard surface agenda (i.e., make use of surface_rtprop_agenda)
ws.Copy(ws.iy_surface_agenda, ws.iy_surface_agenda__UseSurfaceRtprop)
# sensor-only path
ws.Copy(ws.ppath_agenda, ws.ppath_agenda__FollowSensorLosPath)
# no refraction
ws.Copy(ws.ppath_step_agenda, ws.ppath_step_agenda__GeometricPath)
# (standard) emission calculation
ws.Copy(ws.iy_main_agenda, ws.iy_main_agenda__Emission)
# No jacobian calculation
ws.jacobianOff()
# Clearsky = No scattering
ws.cloudboxOff()
# A dummy, empty atmosphere
#
ws.abs_speciesSet(species=[])


@arts_agenda
def propmat_clearsky_agenda(ws):
    ws.Ignore(ws.rtp_mag)
    ws.Ignore(ws.rtp_los)
    ws.Ignore(ws.rtp_pressure)
    ws.Ignore(ws.rtp_temperature)
    ws.Ignore(ws.rtp_vmr)
    ws.Ignore(ws.rtp_nlte)
    ws.Ignore(ws.jacobian_quantities)
    ws.Touch(ws.nlte_source)
    ws.Touch(ws.dpropmat_clearsky_dx)
    ws.Touch(ws.dnlte_dx_source)
    ws.Touch(ws.nlte_dsource_dx)
    ws.propmat_clearskyZero()


ws.propmat_clearsky_agenda = propmat_clearsky_agenda

#
ws.AtmosphereSet1D()
ws.VectorNLogSpace(ws.p_grid, 2, 101300.0, 10000.0)
ws.AtmRawRead(basename="testdata/tropical")
ws.AtmFieldsCalc()
ws.Extract(ws.z_surface, ws.z_field, 0)
ws.Extract(ws.t_surface, ws.t_field, 0)
# Number of Stokes components to be computed
#
ws.IndexSet(ws.stokes_dim, 4)
# A single frequency
#
ws.VectorSet(ws.f_grid, np.array([3.0e10]))
# Set surface to be flat water
#
ws.VectorCreate("data_f_grid")
ws.VectorCreate("data_T_grid")
ws.VectorNLinSpace(ws.data_f_grid, 10, 10000000000.0, 100000000000.0)
ws.VectorNLinSpace(ws.data_T_grid, 5, 270.0, 310.0)
ws.complex_refr_indexWaterLiebe93(
    ws.surface_complex_refr_index, ws.data_f_grid, ws.data_T_grid
)


@arts_agenda
def surface_rtprop_agenda(ws):
    ws.specular_losCalc()
    ws.InterpSurfaceFieldToPosition(out=ws.surface_skin_t, field=ws.t_surface)
    ws.surfaceFlatRefractiveIndex()


ws.surface_rtprop_agenda = surface_rtprop_agenda

# Rayleigh-Jeans brightness temperatures
#
ws.StringSet(ws.iy_unit, "RJBT")
# Create a small across-track scanning sequence as one mblock
#
# We take nadir as reference for our angles.
# It is allowed that sensor_los+mblock_grid is > 180, this is automatically
# mapped to a correct angle.
#
ws.MatrixSet(ws.sensor_pos, np.array([[800000.0]]))
ws.MatrixSet(ws.sensor_los, np.array([[180.0]]))
#
ws.VectorCreate("angles")
ws.VectorNLinSpace(ws.angles, 5, 0.0, 60.0)
ws.IndexCreate("nang")
ws.nelemGet(ws.nang, ws.angles)
#
ws.IndexSet(ws.sensor_norm, 1)
# Here just a dummy value
ws.AntennaOff()
# Simplest way to set some dummy variables
#
ws.Matrix1ColFromVector(ws.mblock_dlos_grid, ws.angles)
# Perform basic checks
ws.atmfields_checkedCalc()
ws.atmgeom_checkedCalc()
ws.cloudbox_checkedCalc()
# Perform RT calculations, with no rotation
#
ws.sensor_responseInit()
ws.sensor_checkedCalc()
ws.yCalc()
# Print(y)
# Perform RT calculations, with rotation following scan angle
# (rotation angle here set to 50% of angle from nadir!)
ws.Copy(ws.stokes_rotation, ws.angles)
ws.VectorScale(ws.stokes_rotation, ws.stokes_rotation, 0.5)
#
ws.sensor_responseInit()
ws.sensor_responseStokesRotation()
ws.sensor_checkedCalc()
ws.yCalc()
# Print(y)
# Check last results
ws.VectorCreate("yref")
ws.VectorSet(
    ws.yref,
    np.array(
        [
            127.654,
            0.0,
            0.0,
            0.0,
            127.694,
            3.97,
            -1.056,
            0.0,
            128.398,
            15.195,
            -8.77282,
            0.0,
            132.644,
            32.479,
            -32.479,
            0.0,
            154.624,
            58.901,
            -102.019,
            0.0,
        ]
    ),
)
ws.Compare(ws.y, ws.yref, 0.05)
