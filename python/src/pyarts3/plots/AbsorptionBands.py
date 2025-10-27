from pyarts3.arts import AbsorptionBands, AscendingGrid, AtmPoint, AtmField, Propmat, SpeciesEnum, PropagationPathPoint
import pyarts3 as pyarts
import matplotlib
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(absorption_bands: AbsorptionBands,
         *,
         mode='normal',
         fig=None,
         ax=None,
         freqs: AscendingGrid | int = 1000,
         atm: AtmPoint | AtmField | None = None,
         pol: Propmat = Propmat(1),
         species: SpeciesEnum = SpeciesEnum.AIR,
         path_point: PropagationPathPoint = PropagationPathPoint(),
         min_pm: float = None,
         cm: matplotlib.colors.ListedColormap = matplotlib.colormaps.get_cmap(
             'viridis'),
         **kwargs):
    """Creates a plot of the absorption by the bands in the AbsorptionBands object.

    If freqs is and index, this is the number of frequency points between the lowest
    and highest frequency in the AbsorptionBands object.  If freqs is an AscendingGrid,
    it is used directly for the x-axis.

    If atm is None, and atmospheric point is created by reading the tropical standard atmosphere
    from the ARTS database and extracting the AtmPoint at position by the path_point.  If atm is
    an AtmField, the AtmPoint at position by the path_point is extracted.

    pol is the dot product of the generated propagation matrix and defaults to the purely
    unpolarized case.

    The species is passed directly to the AbsorptionBands.propagation_matrix() method and defaults
    to all species (AIR).

    The path_point is passed directly to the AbsorptionBands.propagation_matrix() method and defaults
    to pos: [0, 0, 0], los: [0, 0].

    The mode parameter controls the style of the plot.  In 'normal' mode, all absorption lines are plotted
    as the sum of the call to AbsorptionBands.propagation_matrix().  In 'important Y X' modes, the full
    absorption is still shown, but the absorption bands object are split by X, where X is one of 'bands',
    'species', or 'isotopes'.  Y is either 'fill' or nothing.  If it is 'fill', matplotlib's
    fill_between() is used to color the regions, otherwise line plots are used to show the absorption
    of the important contributors.

    Parameters
    ----------
    absorption_bands : ~pyarts3.arts.AbsorptionBands
        The AbsorptionBands object containing the data to plot.
    mode : str, optional
        The mode to use for the plot - see text for descriptions. Defaults to 'normal'.
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    freqs : ~pyarts3.arts.AscendingGrid | int, optional
        The frequency grid to use for the x-axis. Defaults to None manually creating a range.
    atm : ~pyarts3.arts.AtmPoint | ~pyarts3.arts.AtmField | None, optional
        The atmospheric point data to use for the plot. Defaults to None computed values.
    pol : ~pyarts3.arts.Propmat, optional
        The polarization state to consider. Defaults to [1, 0, 0, 0, 0, 0, 0].
    species : ~pyarts3.arts.SpeciesEnum, optional
        The species to include in the plot. Defaults to AIR for all species.
    path_point : ~pyarts3.arts.PropagationPathPoint, optional
        The propagation path point to use for the plot. Defaults to pos: [0, 0, 0], los: [0, 0].
    min_pm : float, optional
        The minimum absorption to consider a band important.  Only used in 'important Y X' modes.
        Defaults to 1% of the minimum absorption in the full propagation matrix.
    cm : ~matplotlib.colors.ListedColormap, optional
        The matplotlib colormap to use for coloring in "important Y X" mode  Must be resampable.  Defaults to 'viridis'.
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig : As input
        As input.
    ax : As input
        As input.
    """
    if isinstance(freqs, int):
        f0 = 1e99
        f1 = -1e99
        for band in absorption_bands:
            for line in absorption_bands[band].lines:
                f0 = min(f0, line.f0)
                f1 = max(f1, line.f0)
        freqs = AscendingGrid(np.linspace(f0, f1, freqs))

    if atm is None:
        ws = pyarts.Workspace()
        ws.absorption_speciesSet(species=[f"{band.isot}" for band in absorption_bands])
        basename = "planets/Earth/afgl/tropical/"
        toa = 1 + path_point.pos[0]
        ws.atmospheric_fieldRead(toa=toa, basename=basename, missing_is_zero=1)
        atm = ws.atmospheric_field(*path_point.pos)
    elif isinstance(atm, AtmField):
        atm = atm(*path_point.pos)

    fig, ax = default_fig_ax(fig, ax, 1, 1, fig_kwargs={'figsize': (6, 4)})

    pm = absorption_bands.propagation_matrix(
        f=freqs, atm=atm, spec=species, path_point=path_point)
    pm = np.einsum('ij,j->i', pm, pol)

    # Plots a simple line plot of the total absorption
    if mode == 'normal':
        select_flat_ax(ax, 0).plot(freqs, pm, color='black', **kwargs)

    # Separates the absorption into important bands/species/isotopes
    elif mode.startswith('important'):
        pm_bands = {}
        keys = None

        # Extracts keys based on mode, computing their propagation matrices on-the-fly
        if mode.endswith('isotopes') or \
           mode.endswith('species'):
            keys = [band.isot if mode.endswith(
                'isotopes') else band.isot.spec for band in absorption_bands]
            keys = list(set(keys))
            for key in keys:
                bands = absorption_bands.extract_species(key)
                pm_band = bands.propagation_matrix(
                    f=freqs, atm=atm, spec=species, path_point=path_point)
                pm_bands[key] = np.einsum('ij,j->i', pm_band, pol)
        elif mode.endswith('bands'):
            keys = [band for band in absorption_bands]
            for key in keys:
                bands = AbsorptionBands({key: absorption_bands[key]})
                pm_band = bands.propagation_matrix(
                    f=freqs, atm=atm, spec=species, path_point=path_point)
                pm_bands[key] = np.einsum('ij,j->i', pm_band, pol)

        min_pm = 0.01 * min(pm) if min_pm is None else min_pm

        # Filter out bands below the minimum propagation matrix value (and set values to NaN)
        for key in keys:
            if max(pm_bands[key]) < min_pm:
                del pm_bands[key]
            else:
                pm_bands[key][pm_bands[key] < min_pm] = np.nan

        # Sort bands so stronger are first (this helps with visibility when plotting)
        pm_keys = list(pm_bands.keys())
        order = []
        for key in pm_keys:
            order.append(np.nansum(pm_bands[key]))
        order = np.array(order)
        order = np.argsort(order)[::-1]
        keys = np.array(pm_keys)[order]

        colors = cm.resampled(len(pm_bands))

        if "fill" in mode:
            select_flat_ax(ax, 0).fill_between(
                freqs, 0, pm, color='black', label='Total', **kwargs)
            for i, key in enumerate(keys):
                select_flat_ax(ax, 0).fill_between(freqs, 0, pm_bands[key],
                                                   color=colors.colors[i], label=str(key), **kwargs)
        else:
            select_flat_ax(ax, 0).plot(
                freqs, pm, color='black', label='Total', **kwargs)
            for i, key in enumerate(keys):
                select_flat_ax(ax, 0).plot(freqs, pm_bands[key],
                                           color=colors.colors[i], label=str(key), **kwargs)

    return fig, ax
