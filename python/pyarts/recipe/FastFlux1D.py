import pyarts

import numpy as np


class FastFlux1D:
    """Creates a Disort clearsky flux operator using the Czarnecki-scheme."""

    def __init__(
        self,
        visibile_surface_reflectivity: float = 0.3,
        thermal_surface_reflectivity: float = 0.05,
        atmospheric_altitude: float = 50e3,
        surface_temperature: float = 300.0,
        max_level_step: float = 1e3,
        NQuad: int = 16,
        atm_latitude: float = 0.0,
        atm_longitude: float = 0.0,
        solar_latitude: float = 0.0,
        solar_longitude: float = 0.0,
    ):
        """Initialization

        Parameters
        ----------
        """

        self.visibile_surface_reflectivity = visibile_surface_reflectivity
        self.thermal_surface_reflectivity = thermal_surface_reflectivity

        self.ws = pyarts.Workspace()

        self.ws.disort_quadrature_dimension = NQuad
        self.ws.disort_fourier_mode_dimension = 1
        self.ws.disort_legendre_polynomial_dimension = 1

        self.ws.absorption_speciesSet(
            species=[
                "H2O-161, H2O-SelfContCKDMT350, H2O-ForeignContCKDMT350",
                "O2-66, O2-CIAfunCKDMT100",
                "N2-44, N2-CIAfunCKDMT252, N2-CIArotCKDMT252",
                "CO2-626, CO2-CKDMT252",
                # "CH4",
                # "O3",
                "O3-XFIT",
            ]
        )

        self.ws.ReadCatalogData()

        for band in self.ws.absorption_bands:
            band.data.cutoff = "ByLine"
            band.data.cutoff_value = 750e9

        self.ws.propagation_matrix_agendaAuto()

        self.ws.surface_fieldSetPlanetEllipsoid(option="Earth")
        self.ws.surface_field["t"] = surface_temperature
        self.ws.absorption_bandsSelectFrequency(fmin=40e9, by_line=1)

        self.ws.atmospheric_fieldRead(
            toa=atmospheric_altitude,
            basename="planets/Earth/afgl/tropical/",
            missing_is_zero=1,
        )

        self.ws.ray_pathGeometricUplooking(
            latitude=atm_latitude,
            longitude=atm_longitude,
            max_step=max_level_step,
        )

        self.visf = pyarts.arts.AscendingGrid.fromxml(
            "planets/Earth/Optimized-Flux-Frequencies/SW-flux-optimized-f_grid.xml"
        )

        self.ir_f = pyarts.arts.AscendingGrid.fromxml(
            "planets/Earth/Optimized-Flux-Frequencies/LW-flux-optimized-f_grid.xml"
        )

        self.visw = pyarts.arts.Vector.fromxml(
            "planets/Earth/Optimized-Flux-Frequencies/SW-flux-optimized-quadrature_weights.xml"
        )
        self.ir_w = pyarts.arts.Vector.fromxml(
            "planets/Earth/Optimized-Flux-Frequencies/LW-flux-optimized-quadrature_weights.xml"
        )

        tmp = pyarts.arts.GriddedField2.fromxml(
            "star/Sun/solar_spectrum_QUIET.xml"
        )
        
        self.ws.sunFromGrid(
            frequency_grid=self.visf,
            sun_spectrum_raw=tmp,
            latitude=solar_latitude,
            longitude=solar_longitude,
        )

    def __call__(
        self,
        atmospheric_profile: dict,
        surface_temperature: float = None,
    ):
        """Call operator to return a propagation matrix

        Parameters
        ----------
        """

        if surface_temperature is not None:
            self.ws.surface_field["t"] = surface_temperature

        N = len(self.ws.ray_path_atmospheric_point)
        for key in atmospheric_profile:
            data = np.atleast_1d(atmospheric_profile[key]).flatten()

            if len(data) == 1:
                for i in range(N):
                    self.ws.ray_path_atmospheric_point[i][key] = data[0]
            elif len(data) == N:
                for i in range(N):
                    self.ws.ray_path_atmospheric_point[i][key] = data[i]
            else:
                assert (
                    False
                ), f"Key: {key}.  Neither size {N} nor 1.  Size: {len(data)}"

        # Visible
        self.ws.frequency_grid = self.visf
        self.ws.ray_path_frequency_gridFromPath()
        self.ws.ray_path_propagation_matrixFromPath()
        self.ws.ray_path_frequency_gridFromPath()
        self.ws.ray_path_propagation_matrixFromPath()
        self.ws.disort_settingsInit()
        self.ws.disort_settingsOpticalThicknessFromPath()
        self.ws.disort_settingsNoLayerThermalEmission()
        self.ws.disort_settingsSurfaceEmissionByTemperature(
            ray_path_point=self.ws.ray_path[0], in_radiance=0
        )
        self.ws.disort_settingsNoSpaceEmission()
        self.ws.disort_settingsSurfaceLambertian(
            value=self.visibile_surface_reflectivity
        )
        self.ws.disort_settingsNoSingleScatteringAlbedo()
        self.ws.disort_settingsNoFractionalScattering()
        self.ws.disort_settingsNoLegendre()
        self.ws.disort_settingsSetSun(
            ray_path_point=self.ws.ray_path[-1], in_radiance=0
        )
        self.ws.disort_spectral_flux_fieldCalc()

        self.SOLAR = np.einsum(
            "i,ijk->jk", self.visw, self.ws.disort_spectral_flux_field
        )

        # IR
        self.ws.frequency_grid = self.ir_f
        self.ws.ray_path_frequency_gridFromPath()
        self.ws.ray_path_propagation_matrixFromPath()
        self.ws.ray_path_frequency_gridFromPath()
        self.ws.ray_path_propagation_matrixFromPath()
        self.ws.disort_settingsInit()
        self.ws.disort_settingsOpticalThicknessFromPath()
        self.ws.disort_settingsLayerThermalEmissionLinearInTau(in_radiance=0)
        self.ws.disort_settingsSurfaceEmissionByTemperature(
            ray_path_point=self.ws.ray_path[0], in_radiance=0
        )
        self.ws.disort_settingsCosmicMicrowaveBackgroundRadiation(
            in_radiance=0
        )
        self.ws.disort_settingsSurfaceLambertian(
            value=self.thermal_surface_reflectivity
        )
        self.ws.disort_settingsNoSingleScatteringAlbedo()
        self.ws.disort_settingsNoFractionalScattering()
        self.ws.disort_settingsNoLegendre()
        self.ws.disort_settingsNoSun()
        self.ws.disort_spectral_flux_fieldCalc()

        self.THERMAL = np.einsum(
            "i,ijk->jk", self.ir_w, self.ws.disort_spectral_flux_field
        )

        return 0.0
