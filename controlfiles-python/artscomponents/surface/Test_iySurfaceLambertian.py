#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Manfred Brath
"""

import numpy as np
from pyarts.workspace import Workspace, arts_agenda
from pyarts.arts import constant
from pyarts import xml


# =============================================================================
# %%  define scattering agenda
# =============================================================================


# surface scattering agenda
@arts_agenda
def iy_surface_agenda(ws):
    ws.iySurfaceInit()
    ws.Ignore(ws.dsurface_rmatrix_dx)
    ws.Ignore(ws.dsurface_emission_dx)

    ws.iySurfaceLambertian()
    ws.iySurfaceLambertianDirect()


# =============================================================================
# %%  define RT functions
# =============================================================================


def ARTS_clearsky(f_grid, sensor_pos, sensor_los, sun_longitude_pos,
                  surface_scalar_reflectivity, npres=81, nlat=3, nlon=5, stokes_dim=1):
    '''
    Calculates the clear sky spectral radiance. As no scattering is present,
    radiation is seen when looking directly towards the sun or when looking
    towards its specular reflection from the ground.

    Args:
        f_grid (ndarray): Frequncy grid in Hz.
        sensor_pos (ndarray): Sensor position [Altitude, Latitude, Longitude].
        sensor_los (ndarray): Sensor looking direction [zenith angle, azimuth angle].
        sun_longitude_pos (float): Zenith longitude position of the sun.
        ComplexRefractiveIndex (ndarray): Surface complex refractive index.
        npres (int, optional): Number of pressure levels. Defaults to 81.
        nlat (int, optional): Number of latitude points. Defaults to 3.
        nlon (int, optional): Number of longitude points. Defaults to 5.
        stokes_dim (int, optional): Stoke dimension. Defaults to 2.

    Returns:
        ndarray: Spectral radiance at sensor.

    '''

    # =============================================================================
    # open workspace
    # =============================================================================

    ws = Workspace()
    ws.verbositySetScreen(level=2)

    # import basic definitions
    ws.execute_controlfile("general/continua.arts")
    ws.execute_controlfile("general/planet_earth.arts")

    # =============================================================================
    # select/define agendas
    # =============================================================================

    # set agenda
    ws.iy_main_agendaSet( option="Clearsky" )
    ws.iy_space_agendaSet( option="CosmicBackground" )
    ws.ppath_agendaSet( option="FollowSensorLosPath" )
    ws.ppath_step_agendaSet( option="GeometricPath" )
    ws.water_p_eq_agendaSet()

    ws.iy_surface_agenda = iy_surface_agenda

    # =============================================================================
    # basic conditions
    # =============================================================================

    # Postion and line-of-sight of sensor
    ws.sensor_pos = sensor_pos
    ws.sensor_los = sensor_los
    ws.VectorSet(ws.rte_pos2, [])

    # =============================================================================
    # define environment
    # =============================================================================

    # Number of Stokes components to be computed
    #
    ws.IndexSet(ws.stokes_dim, stokes_dim)

    # Reference ellipsoid
    ws.refellipsoidEarth(ws.refellipsoid, "Sphere")

    # Frequency grid
    ws.VectorSet(ws.f_grid, f_grid)

    # A pressure grid rougly matching 0 to 80 km, in steps of 2 km.
    ws.VectorNLogSpace(ws.p_grid, npres, 1013e2, 1)

    # Atmospheric dimensionality and lat/lon grids
    ws.VectorNLinSpace(ws.lat_grid, nlat, -90., 90.)
    ws.VectorNLinSpace(ws.lon_grid, nlon, -180., 180.)
    ws.AtmosphereSet3D()

    # Definition of species
    ws.abs_speciesSet(species=
                      ["N2",
                       "O2"])

    # No line data needed here
    ws.abs_lines_per_speciesSetEmpty()

    ws.propmat_clearsky_agendaAuto()

    # Atmospheric profiles
    ws.AtmRawRead(basename="testdata/tropical")
    ws.AtmFieldsCalcExpand1D()

    # Get ground altitude (z_surface) from z_field
    ws.MatrixSetConstant(ws.z_surface, nlat, nlon, 0.)

    # set surface temperature
    ws.MatrixSetConstant(ws.t_surface, nlat, nlon, 1.)
    # ws.surface_skin_t=ws.t_surface.value[0,0]

    ws.ArrayOfStringSet(ws.surface_props_names, ["Skin temperature"])
    ws.Tensor3SetConstant(ws.surface_props_data, 1, nlat, nlon, ws.t_surface.value[0, 0])

    ws.surface_scalar_reflectivity = surface_scalar_reflectivity

    # set star source
    ws.starsAddSingleBlackbody(longitude=sun_longitude_pos)

    # =============================================================================
    # the calculation
    # =============================================================================

    # No jacobian calculations
    ws.jacobianOff()

    # No particulate scattering
    ws.cloudboxOff()

    # No gas scattering
    ws.gas_scatteringOff()

    # No sensor
    ws.sensorOff()

    # Check model atmosphere
    ws.atmfields_checkedCalc()
    ws.atmgeom_checkedCalc()
    ws.cloudbox_checkedCalc()
    ws.sensor_checkedCalc()

    # Propagation path agendas and variables
    ws.NumericSet(ws.ppath_lmax, -1)

    ws.lbl_checkedCalc()

    # Switch off stars
    ws.IndexSet(ws.stars_do, 1)

    ws.yCalc()

    y = ws.y.value.value

    return y, ws


# =============================================================================
# %% paths/constants
# =============================================================================
if __name__ == "__main__":
    # Sensor position in altitude/lat/lon
    sensor_pos = np.array([[2., 0., 0.]])

    # Wavelength im nm
    wavelength = np.array([600]) * 1e-9

    # Frequency grid
    f_grid = constant.c / wavelength

    # Set Complex refractive index similar to water at 600 nm
    surface_reflectivity = np.array([1.])

    # Set the sun zenith longitude position. Here we only need the the longitude,
    # because by default sun's zenith is by default at 0 deg latitude.
    sun_longitude_pos = 45.

    # Sensor looking direction in zenith angle (0 = upwards, 180 = downward) and
    # azimuth angle ( 0 = North, 90 = east)
    sensor_los = np.array([[180. - sun_longitude_pos, 90.]])

    y, ws = ARTS_clearsky(f_grid, sensor_pos, sensor_los, sun_longitude_pos,
                          surface_reflectivity)

    # Reference data
    yREFERENCE = np.array([4.49321403187969e-13])

    # Compare with reference
    ws.CompareRelative(y, yREFERENCE, 1e-6)
