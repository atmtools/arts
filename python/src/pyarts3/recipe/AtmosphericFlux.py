import pyarts3 as pyarts

import numpy as np

from dataclasses import dataclass


@dataclass(slots=True)
class Flux:
    name: str
    up: np.ndarray
    diffuse_down: np.ndarray
    direct_down: np.ndarray

    @property
    def down(self):
        return self.diffuse_down + self.direct_down


class AtmosphericFlux:
    """Creates a Disort clearsky flux operator using the Czarnecki-scheme."""

    def __init__(
        self,
        visible_surface_reflectivity: float = 0.3,
        thermal_surface_reflectivity: float = 0.05,
        atmospheric_altitude: float = 50e3,
        surface_temperature: float = 300.0,
        max_level_step: float = 1e3,
        NQuad: int = 16,
        atm_latitude: float = 0.0,
        atm_longitude: float = 0.0,
        solar_latitude: float = 0.0,
        solar_longitude: float = 0.0,
        species: list = ["H2O-161", "O2-66", "N2-44", "CO2-626", "O3-XFIT"],
        remove_lines_percentile: dict[pyarts.arts.SpeciesEnum,
                                      float] | float | None = None,
    ):
        """Compute the total flux for a given atmospheric profile and surface temperature

        The operator allows you to change the

        Args:
            visible_surface_reflectivity (float, optional): The surface reflectivity constant for Disort in visible. Defaults to 0.3.
            thermal_surface_reflectivity (float, optional): The surface reflectivity constant for Disort in thermal. Defaults to 0.05.
            atmospheric_altitude (float, optional): The top-of-the-atmosphere altitude [m]. Defaults to 50e3.
            surface_temperature (float, optional): The surface temperature [K]. Defaults to 300.0.
            max_level_step (float, optional): The maximum thickness of layers [m]. Defaults to 1e3.
            NQuad (int, optional): The number of quadratures used by Disort. Defaults to 16.
            atm_latitude (float, optional): Latitude of profile [degrees]. Defaults to 0.0.
            atm_longitude (float, optional): Longitude of profile [degrees]. Defaults to 0.0.
            solar_latitude (float, optional): Latitude of sun [degrees]. Defaults to 0.0.
            solar_longitude (float, optional): Longitude of sun [degrees]. Defaults to 0.0.
            species (list, optional): The list of absorption species. Defaults to [ "H2O-161", "O2-66", "N2-44", "CO2-626", "O3-XFIT", ].
            remove_lines_percentile (dict | float | None, optional): The percentile of lines to remove [0, 100].  Per species if dict. Defaults to None.
        """

        self.visible_surface_reflectivity = visible_surface_reflectivity
        self.thermal_surface_reflectivity = thermal_surface_reflectivity

        self.ws = pyarts.Workspace()

        self.ws.disort_quadrature_dimension = NQuad
        self.ws.disort_fourier_mode_dimension = 1
        self.ws.disort_legendre_polynomial_dimension = 1

        self.ws.absorption_speciesSet(species=species)

        self.ws.ReadCatalogData()

        for band in self.ws.absorption_bands:
            self.ws.absorption_bands[band].cutoff = "ByLine"
            self.ws.absorption_bands[band].cutoff_value = 750e9

        if remove_lines_percentile is not None:
            self.ws.absorption_bands.keep_hitran_s(remove_lines_percentile)

        self.ws.propagation_matrix_agendaAuto()

        self.ws.surface_fieldPlanet(option="Earth")
        self.ws.surface_field["t"] = surface_temperature
        self.ws.absorption_bandsSelectFrequencyByLine(fmin=40e9)

        self.ws.atmospheric_fieldRead(
            toa=atmospheric_altitude,
            basename="planets/Earth/afgl/tropical/",
            missing_is_zero=1,
        )

        self.ws.ray_pathGeometricDownlooking(
            latitude=atm_latitude,
            longitude=atm_longitude,
            max_stepsize=max_level_step,
        )

        self.ws.ray_path_atmospheric_pointFromPath()

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

        tmp = pyarts.arts.GriddedField2.fromxml("star/Sun/solar_spectrum_QUIET.xml")
        self.ws.sunFromGrid(
            frequency_grid=self.visf,
            sun_spectrum_raw=tmp,
            latitude=solar_latitude,
            longitude=solar_longitude,
        )

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
            self.ws.ray_path_atmospheric_point.to_dict(
                core=core, specs=specs, nlte=nlte, ssprops=ssprops, isots=isots
            )
        )

    def __call__(
        self,
        atmospheric_profile: dict = {},
        surface_temperature: float = None,
    ):
        """Get the total flux profile

        Args:
            atmospheric_profile (dict, optional): A dictionary of atmospheric data. Defaults to {}.
            surface_temperature (float, optional): A surface temperature. Defaults to None.

        Returns:
            Flux, Flux, numpy.ndarray: The solar and thermal fluxes and the center altitudes of the layers.
        """

        if surface_temperature is not None:
            self.ws.surface_field["t"] = surface_temperature

        self.ws.ray_path_atmospheric_point.update(atmospheric_profile)

        # Visible
        self.ws.frequency_grid = self.visf
        self.ws.disort_settings_agendaSetup(
            layer_emission_setting="None",
            scattering_setting="None",
            space_setting="None",
            surface_setting="Lambertian",
            sun_setting="Sun",
            surface_lambertian_value=self.visible_surface_reflectivity *
            np.ones_like(self.visf)
        )
        self.ws.disort_spectral_flux_fieldFromAgenda()

        self.SOLAR = Flux("solar",
                          np.einsum("i,ik->k", self.visw,
                                    self.ws.disort_spectral_flux_field.up),
                          np.einsum("i,ik->k", self.visw,
                                    self.ws.disort_spectral_flux_field.down_diffuse),
                          np.einsum("i,ik->k", self.visw, self.ws.disort_spectral_flux_field.down_direct))

        # IR
        self.ws.frequency_grid = self.ir_f
        self.ws.disort_settings_agendaSetup(
            layer_emission_setting="LinearInTau",
            scattering_setting="None",
            space_setting="CosmicMicrowaveBackgroundRadiation",
            surface_setting="ThermalLambertian",
            sun_setting="None",
            surface_lambertian_value=self.thermal_surface_reflectivity *
            np.ones_like(self.visf)
        )
        self.ws.disort_spectral_flux_fieldFromAgenda()

        self.THERMAL = Flux("thermal",
                            np.einsum("i,ik->k", self.ir_w,
                                      self.ws.disort_spectral_flux_field.up),
                            np.einsum("i,ik->k", self.ir_w,
                                      self.ws.disort_spectral_flux_field.down_diffuse),
                            np.einsum("i,ik->k", self.ir_w, self.ws.disort_spectral_flux_field.down_direct))

        return (
            self.SOLAR,
            self.THERMAL,
            np.array(
                [
                    0.5 * (self.ws.ray_path[i].pos[0] + self.ws.ray_path[i + 1].pos[0])
                    for i in range(len(self.ws.ray_path) - 1)
                ]
            ),
        )
