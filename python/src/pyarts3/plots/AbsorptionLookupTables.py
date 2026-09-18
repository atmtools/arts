"""Plot absorption lookup-table cross-sections and column opacities."""

from collections.abc import Mapping, Sequence
from math import ceil

import matplotlib
import numpy as np
from pyarts3.arts import AbsorptionLookupTables, SpeciesEnum, constants


from .common import default_fig_ax, select_flat_ax

__all__ = ["plot"]

_BOLTZMANN = constants.k
_DRY_AIR_GAS_CONSTANT = 287.058
_STANDARD_GRAVITY = 9.80665


def _selected_species(
    data: AbsorptionLookupTables,
    species: SpeciesEnum | str | Sequence[SpeciesEnum | str] | None,
) -> list[SpeciesEnum]:
    available = list(data)
    if species is None:
        return available
    if isinstance(species, (SpeciesEnum, str)):
        species = [species]

    selected = []
    by_name = {str(key): key for key in available}
    for requested in species:
        name = str(requested)
        if name not in by_name:
            raise KeyError(
                f"Species {name!r} is not in the lookup table. "
                f"Available species: {', '.join(by_name)}"
            )
        selected.append(by_name[name])
    return selected


def _perturbation_index(grid, index: int, name: str) -> int:
    if grid is None or not len(grid):
        return 0
    size = len(grid)
    if not -size <= index < size:
        raise IndexError(f"{name} index {index} is outside a grid of size {size}")
    return index % size


def _table_slice(table, temperature_perturbation: int, water_perturbation: int):
    it = _perturbation_index(
        table.t_pert, temperature_perturbation, "Temperature perturbation"
    )
    iw = _perturbation_index(table.w_pert, water_perturbation, "Water perturbation")
    return np.asarray(table.xsec)[it, iw], it, iw


def _profile_for_species(profiles: Mapping | None, species: SpeciesEnum):
    if profiles is None:
        return None
    if species in profiles:
        return profiles[species]
    if str(species) in profiles:
        return profiles[str(species)]
    return None


def _altitude_from_pressure(pressure, temperature, gravity):
    layer_temperature = 0.5 * (temperature[:-1] + temperature[1:])
    layer_height = (
        _DRY_AIR_GAS_CONSTANT
        * layer_temperature
        / gravity
        * np.log(pressure[:-1] / pressure[1:])
    )
    return np.concatenate(([0.0], np.cumsum(layer_height)))


def _opacity(
    data,
    selected,
    temperature_perturbation,
    water_perturbation,
    vmr_profiles,
    altitude,
    gravity,
):
    result = {}
    common_frequency = None

    for species in selected:
        table = data[species]
        xsec, it, iw = _table_slice(table, temperature_perturbation, water_perturbation)
        frequency = np.asarray(table.f_grid)
        pressure = np.exp(np.asarray(table.log_p_grid))
        temperature = np.asarray(table.t_atmref, dtype=float).copy()
        if table.t_pert is not None and len(table.t_pert):
            temperature += np.asarray(table.t_pert)[it]

        vmr = _profile_for_species(vmr_profiles, species)
        if vmr is None and species == SpeciesEnum.H2O:
            vmr = np.asarray(table.water_atmref, dtype=float)
            if table.w_pert is not None and len(table.w_pert):
                vmr = vmr * np.asarray(table.w_pert)[iw]
        if vmr is None:
            raise ValueError(
                f"Opacity for {species} requires vmr_profiles[{str(species)!r}]"
            )
        vmr = np.broadcast_to(np.asarray(vmr, dtype=float), pressure.shape)

        height = (
            _altitude_from_pressure(pressure, temperature, gravity)
            if altitude is None
            else np.broadcast_to(np.asarray(altitude, dtype=float), pressure.shape)
        )
        number_density = pressure * vmr / (_BOLTZMANN * temperature)
        result[species] = np.trapezoid(xsec * number_density[:, None], height, axis=0)

        if common_frequency is None:
            common_frequency = frequency
        elif not np.array_equal(common_frequency, frequency):
            raise ValueError("Opacity plotting requires a common frequency grid")

    return common_frequency, result


def plot(
    data: AbsorptionLookupTables,
    *,
    mode: str = "cross_section",
    fig: matplotlib.figure.Figure | None = None,
    ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | np.ndarray | None = None,
    species: SpeciesEnum | str | Sequence[SpeciesEnum | str] | None = None,
    pressures: Sequence[float] | None = None,
    temperature_perturbation: int = 0,
    water_perturbation: int = 0,
    vmr_profiles: Mapping[SpeciesEnum | str, Sequence[float] | float] | None = None,
    altitude: Sequence[float] | None = None,
    gravity: float = _STANDARD_GRAVITY,
    cols: int = 3,
    total: bool = True,
    unity_line: bool = True,
    **kwargs,
) -> tuple[
    matplotlib.figure.Figure,
    matplotlib.axes.Axes | list[matplotlib.axes.Axes] | np.ndarray,
]:
    """Plot cross-sections or column opacities from ARTS lookup tables.

    In ``"cross_section"`` mode, one panel is created per species and one line
    per selected pressure.  In ``"opacity"`` mode, the vertically integrated
    opacity of each species is plotted on one panel.  Opacity requires a VMR
    profile for every selected species except H2O, whose stored reference
    profile is used by default.

    Parameters
    ----------
    data : ~pyarts3.arts.AbsorptionLookupTables
        Lookup-table map, keyed by species.
    mode : {"cross_section", "opacity"}, optional
        Quantity to plot.
    fig : matplotlib.figure.Figure, optional
        Existing Matplotlib figure.
    ax : matplotlib.axes.Axes or numpy.ndarray, optional
        Existing Matplotlib axes.
    species : SpeciesEnum, str, or sequence, optional
        Species to include.  The default includes all tables.
    pressures : sequence of float, optional
        Pressure levels [Pa] nearest to those plotted in cross-section mode.
        At most six approximately evenly spaced levels are used by default.
    temperature_perturbation, water_perturbation : int, optional
        Indices along the corresponding lookup-table dimensions.
    vmr_profiles : mapping, optional
        Species VMR profiles used for opacity. Scalars are broadcast over the
        pressure grid.
    altitude : sequence of float, optional
        Altitudes [m] corresponding to the lookup pressure grid. If omitted,
        hydrostatic layer thicknesses are estimated from the stored reference
        temperatures.
    gravity : float, optional
        Gravity [m/s2] used for the hydrostatic altitude estimate.
    cols : int, optional
        Maximum panel columns in cross-section mode.
    total, unity_line : bool, optional
        Add total-opacity and opacity-one reference lines in opacity mode.
    **kwargs
        Forwarded to Matplotlib's ``plot`` calls.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The Matplotlib figure.
    ax : matplotlib.axes.Axes or numpy.ndarray
        The Matplotlib axes.
    """
    if not data:
        raise ValueError("Cannot plot empty absorption lookup tables")
    selected = _selected_species(data, species)
    if not selected:
        raise ValueError("No species selected")

    normalized_mode = mode.lower().replace("-", "_")
    if normalized_mode in {"cross_section", "cross_sections", "xsec"}:
        if cols < 1:
            raise ValueError("cols must be positive")
        ncols = min(cols, len(selected))
        nrows = ceil(len(selected) / ncols)
        fig, ax = default_fig_ax(
            fig,
            ax,
            nrows,
            ncols,
            N=len(selected),
            fig_kwargs={"figsize": (4 * ncols, 3 * nrows)},
        )

        for index, selected_species in enumerate(selected):
            table = data[selected_species]
            xsec, _, _ = _table_slice(
                table, temperature_perturbation, water_perturbation
            )
            pressure = np.exp(np.asarray(table.log_p_grid))
            if pressures is None:
                pressure_indices = np.unique(
                    np.linspace(0, len(pressure) - 1, min(6, len(pressure)), dtype=int)
                )
            else:
                pressure_indices = np.asarray(
                    [np.abs(pressure - value).argmin() for value in pressures]
                )

            current_ax = select_flat_ax(ax, index)
            for pressure_index in pressure_indices:
                current_ax.plot(
                    table.f_grid,
                    xsec[pressure_index],
                    label=f"{pressure[pressure_index] / 100:g} hPa",
                    **kwargs,
                )
            current_ax.set_yscale("log")
            current_ax.set_title(str(selected_species))
            current_ax.set_xlabel("Frequency [Hz]")
            current_ax.set_ylabel("Cross-section [m$^2$]")
            current_ax.legend(fontsize="small")
        return fig, ax

    if normalized_mode == "opacity":
        fig, ax = default_fig_ax(fig, ax, fig_kwargs={"figsize": (6, 4)})
        current_ax = select_flat_ax(ax, 0)
        frequency, opacities = _opacity(
            data,
            selected,
            temperature_perturbation,
            water_perturbation,
            vmr_profiles,
            altitude,
            gravity,
        )
        for selected_species, values in opacities.items():
            current_ax.plot(frequency, values, label=str(selected_species), **kwargs)
        if total:
            current_ax.plot(
                frequency,
                np.sum(list(opacities.values()), axis=0),
                color="black",
                linewidth=1,
                label="Total",
            )
        if unity_line:
            current_ax.axhline(1.0, color="black", linewidth=1, linestyle="--")
        current_ax.set_yscale("log")
        current_ax.set_xlabel("Frequency [Hz]")
        current_ax.set_ylabel("Column opacity")
        current_ax.legend(fontsize="small")
        return fig, ax

    raise ValueError("mode must be 'cross_section' (or 'xsec') or 'opacity'")
