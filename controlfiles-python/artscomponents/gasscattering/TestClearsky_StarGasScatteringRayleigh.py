#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Manfred Brath
"""

import numpy as np
from pyarts.workspace import Workspace, arts_agenda



# =============================================================================
# %%  define scattering agenda
# =============================================================================

# gas scattering agenda
@arts_agenda
def gas_scattering_agenda(ws):
    ws.Ignore(ws.rtp_vmr)
    ws.gas_scattering_coefAirSimple()
    ws.gas_scattering_matRayleigh()


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


def starARTS_clearsky(f_grid, sensor_pos, sensor_los, sun_pos, Reflectivity,
                      npres=81, nlat=3, nlon=5, stokes_dim=1):
    '''
    Calculates the clear sky spectral radiance using the above defined gas
    scattering (see gas_scattering_agenda) and surface properties
    (see iy_surface_agenda).

    Args:
        f_grid (ndarray): Frequncy grid in Hz.
        sensor_pos (ndarray): Sensor position [Altitude, Latitude, Longitude].
        sensor_los (ndarray): Sensor looking direction [zenith angle, azimuth angle].
        sun_pos (ndarray): Zenith position of the sun.
        Reflectivity (ndarray): Surface reflectivity.
        npres (int, optional): Number of pressure levels. Defaults to 81.
        nlat (int, optional): Number of latitude points. Defaults to 3.
        nlon (int, optional): Number of longitude points. Defaults to 5.
        stokes_dim (int, optional): Stoke dimension. Defaults to 1.

    Returns:
        ndarray: Spectral radiance at sensor.

    '''

    # =============================================================================
    # open workspace
    # =============================================================================

    ws = Workspace()
    ws.verbositySetScreen(level=2)

    # =============================================================================
    # select/define agendas
    # =============================================================================

    ws.execute_controlfile("general/continua.arts")
    ws.execute_controlfile("general/planet_earth.arts")

    # cosmic background radiation
    ws.iy_space_agendaSet( option="CosmicBackground" )

    # sensor-only path
    ws.ppath_agendaSet( option="FollowSensorLosPath" )

    # no refraction
    ws.ppath_step_agendaSet( option="GeometricPath" )

    # main agenda
    ws.iy_main_agendaSet( option="Clearsky")

    # water agenda
    ws.water_p_eq_agendaSet()

    # surface agenda
    ws.iy_surface_agenda = iy_surface_agenda

    # gas scattering agenda
    ws.gas_scattering_agenda = gas_scattering_agenda

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
    ws.IndexSet(ws.stokes_dim, stokes_dim)

    # Reference ellipsoid
    ws.refellipsoidEarth(ws.refellipsoid, "Sphere")

    # Frequency grid
    ws.VectorSet(ws.f_grid, f_grid)

    # Define pressure grid
    ws.VectorNLogSpace(ws.p_grid, npres, 1013e2, 1)

    # Atmospheric dimensionality and lat/lon grids
    ws.VectorNLinSpace(ws.lat_grid, nlat, -90., 90.)
    ws.VectorNLinSpace(ws.lon_grid, nlon, -180., 180.)
    ws.AtmosphereSet3D()

    # Definition of species
    ws.abs_speciesSet(species=
                      ["H2O",
                       "N2",
                       "O2"])

    # No line data needed here, because here we consider only rayleigh scattering
    ws.abs_lines_per_speciesSetEmpty()

    ws.propmat_clearsky_agendaAuto()

    # Load atmospheric data
    ws.AtmRawRead(basename="testdata/tropical")

    # Expand 1D atmospheric profile to 3D
    ws.AtmFieldsCalcExpand1D()

    # Get ground altitude (z_surface) from z_field
    ws.MatrixSetConstant(ws.z_surface, nlat, nlon, 0.)

    # set surface skin temperature
    ws.ArrayOfStringSet(ws.surface_props_names, ["Skin temperature"])
    ws.Tensor3SetConstant(ws.surface_props_data, 1, nlat, nlon, ws.t_field.value[0, 0, 0])

    # Set surface relectivity
    ws.VectorSet(ws.surface_scalar_reflectivity, Reflectivity)

    # set a simple blackbody sun
    ws.starsAddSingleBlackbody(latitude=sun_pos[0], longitude=sun_pos[1])

    # No jacobian calculations
    ws.jacobianOff()

    # No particulate scattering
    ws.cloudboxOff()

    # No sensor model
    ws.sensorOff()

    # Set ropagation path length
    ws.NumericSet(ws.ppath_lmax, -1)
    ws.NumericSet(ws.ppath_lraytrace, 1e4)

    # =============================================================================
    # the calculation
    # =============================================================================

    # Check model atmosphere
    ws.atmfields_checkedCalc()
    ws.atmgeom_checkedCalc()
    ws.cloudbox_checkedCalc()
    ws.sensor_checkedCalc()
    ws.lbl_checkedCalc()

    # Switch on gas scattering
    ws.IndexSet(ws.gas_scattering_do, 1)

    # Switch off stars
    ws.IndexSet(ws.stars_do, 1)

    # the actual simulation
    ws.yCalc()

    return np.copy(ws.y.value), ws


# =============================================================================
# %%
# =============================================================================

if __name__ == "__main__":
    # Set frequency
    f_grid = np.array([6e14])  # corresponds to wavelength of roughly 500 nm

    # Sensor position in altitude/lat/lon
    sensor_pos = np.array([[1000e3, 0., 0.]])

    # Sensor looking direction in zenith angle (0 = upwards, 180 = downward) and
    # azimuth angle ( 0 = North, 90 = east)
    sensor_los = np.array([[180., 0.]])

    # Sun position in lat/lon. The sun postion defines where the sun is at zenith.
    sun_pos = np.array([0., 0.])

    # Reflectivity
    Reflectivity = np.array([1.])

    y, ws = starARTS_clearsky(f_grid, sensor_pos, sensor_los, sun_pos, Reflectivity, npres=81, nlat=3, nlon=5,
                              stokes_dim=1)

    # Reference data
    yREFERENCE = np.array([3.94192129256117e-13])

    # Compare with reference
    ws.CompareRelative(y, yREFERENCE, 1e-6)
