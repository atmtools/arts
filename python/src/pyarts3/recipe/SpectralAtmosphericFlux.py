import pyarts3 as pyarts

import numpy as np

from dataclasses import dataclass


@dataclass(slots=True)
class Flux:
    up: np.ndarray
    diffuse_down: np.ndarray
    direct_down: np.ndarray

    @property
    def down(self):
        return self.diffuse_down + self.direct_down


class SpectralAtmosphericFlux:
    """Creates a Disort clearsky flux operator using the Czarnecki-scheme."""

    def __init__(
        self,
        visible_surf_reflectivity: float = 0.3,
        thermal_surf_reflectivity: float = 0.05,
        atmospheric_altitude: float = 50e3,
        surf_temperature: float = 300.0,
        max_level_step: float = 1e3,
        NQuad: int = 16,
        atm_latitude: float = 0.0,
        atm_longitude: float = 0.0,
        solar_latitude: float = 0.0,
        solar_longitude: float = 0.0,
        species=["H2O-161", "O2-66", "N2-44", "CO2-626", "O3-XFIT"],
        remove_lines_percentile: dict[pyarts.arts.SpeciesEnum,
                                      float] | float | None = None,
    ):
        """Compute the total flux for a given atmospheric profile and surface temperature

        The operator allows you to change the

        Args:
            visible_surf_reflectivity (float, optional): The surface reflectivity constant for Disort in visible. Defaults to 0.3.
            thermal_surf_reflectivity (float, optional): The surface reflectivity constant for Disort in thermal. Defaults to 0.05.
            atmospheric_altitude (float, optional): The top-of-the-atmosphere altitude [m]. Defaults to 50e3.
            surf_temperature (float, optional): The surface temperature [K]. Defaults to 300.0.
            max_level_step (float, optional): The maximum thickness of layers [m]. Defaults to 1e3.
            NQuad (int, optional): The number of quadratures used by Disort. Defaults to 16.
            atm_latitude (float, optional): Latitude of profile [degrees]. Defaults to 0.0.
            atm_longitude (float, optional): Longitude of profile [degrees]. Defaults to 0.0.
            solar_latitude (float, optional): Latitude of sun [degrees]. Defaults to 0.0.
            solar_longitude (float, optional): Longitude of sun [degrees]. Defaults to 0.0.
            species (list, optional): The list of absorption species. Defaults to [ "H2O-161", "O2-66", "N2-44", "CO2-626", "O3-XFIT", ].
            remove_lines_percentile (dict | float | None, optional): The percentile of lines to remove [0, 100].  Per species if dict. Defaults to None.
        """

        self.visible_surf_reflectivity = visible_surf_reflectivity
        self.thermal_surf_reflectivity = thermal_surf_reflectivity

        self.ws = pyarts.Workspace()

        self.ws.disort_quadrature_dimension = NQuad
        self.ws.disort_fourier_mode_dimension = 1
        self.ws.disort_legendre_polynomial_dimension = 1

        self.ws.abs_speciesSet(species=species)

        self.ws.ReadCatalogData()

        for band in self.ws.abs_bands:
            self.ws.abs_bands[band].cutoff = "ByLine"
            self.ws.abs_bands[band].cutoff_value = 750e9

        if remove_lines_percentile is not None:
            self.ws.abs_bands.keep_hitran_s(remove_lines_percentile)

        self.ws.spectral_propmat_agendaAuto()

        self.ws.surf_fieldPlanet(option="Earth")
        self.ws.surf_field["t"] = surf_temperature
        self.ws.abs_bandsSelectFrequencyByLine(fmin=40e9)

        self.ws.atm_fieldRead(
            toa=atmospheric_altitude,
            basename="planets/Earth/afgl/tropical/",
            missing_is_zero=1,
        )

        self.ws.ray_pathGeometricDownlooking(
            lat=atm_latitude,
            lon=atm_longitude,
            max_stepsize=max_level_step,
        )

        self.ws.atm_pathFromPath()

        self.sun = pyarts.arts.GriddedField2.fromxml(
            "star/Sun/solar_spectrum_QUIET.xml"
        )
        self.solar_latitude = solar_latitude
        self.solar_longitude = solar_longitude

    def get_atmosphere(
        self, core=True, specs=True, nlte=False, ssprops=False, isots=False
    ):
        """Return the atmospheric field as a dictionary of python types.

        Args:
            core (bool, optional): See :meth:`ArrayOfAtmPoint.to_dict`. Defaults to True.
            specs (bool, optional): See :meth:`ArrayOfAtmPoint.to_dict`. Defaults to True.
            nlte (bool, optional): See :meth:`ArrayOfAtmPoint.to_dict`. Defaults to False.
            ssprops (bool, optional): See :meth:`ArrayOfAtmPoint.to_dict`. Defaults to False.
            isots (bool, optional): See :meth:`ArrayOfAtmPoint.to_dict`. Defaults to False.

        Returns:
            dict: Atmospheric field dictionary
        """
        return pyarts.arts.stringify_keys(
            self.ws.atm_path.to_dict(
                core=core, specs=specs, nlte=nlte, ssprops=ssprops, isots=isots
            )
        )

    def __call__(
        self,
        freq_grid: pyarts.arts.AscendingGrid,
        atm_profile: dict = {},
        surf_temperature: float = None,
    ):
        """Get the total flux profile

        Args:
            freq_grid (pyarts3.arts.AscendingGrid): The frequency grid
            atm_profile (dict, optional): The atmospheric profile. Defaults to {}.
            surf_temperature (float, optional): The surface temperature. Defaults to None.

        Returns:
            Flux, numpy.ndarray: Flux profile and average layer altitudes
        """

        if surf_temperature is not None:
            self.ws.surf_field["t"] = surf_temperature

        self.ws.atm_path.update(atm_profile)

        # Visible
        self.ws.freq_grid = freq_grid

        self.ws.sunFromGrid(
            sun_spectrum_raw=self.sun,
            lat=self.solar_latitude,
            lon=self.solar_longitude,
        )

        self.ws.freq_grid_pathFromPath()
        self.ws.spectral_propmat_pathFromPath()
        self.ws.disort_settingsInit()
        self.ws.disort_settingsOpticalThicknessFromPath()
        self.ws.disort_settingsLayerThermalEmissionLinearInTau()
        self.ws.disort_settingsSurfaceEmissionByTemperature(
            ray_path_point=self.ws.ray_path[0]
        )
        self.ws.disort_settingsCosmicMicrowaveBackgroundRadiation()
        self.ws.disort_settingsSurfaceLambertian(
            value=self.visible_surf_reflectivity
        )
        self.ws.disort_settingsNoSingleScatteringAlbedo()
        self.ws.disort_settingsNoFractionalScattering()
        self.ws.disort_settingsNoLegendre()
        self.ws.disort_settingsSetSun(ray_path_point=self.ws.ray_path[-1])
        self.ws.disort_spectral_flux_fieldCalc()

        # Shape is f x 3 x np, we want 3 x f x np

        return (
            Flux(self.ws.disort_spectral_flux_field.up, self.ws.disort_spectral_flux_field.down_diffuse,
                 self.ws.disort_spectral_flux_field.down_direct),
            np.array(
                [
                    0.5 * (self.ws.ray_path[i].pos[0] + self.ws.ray_path[i + 1].pos[0])
                    for i in range(len(self.ws.ray_path) - 1)
                ]
            ),
        )
