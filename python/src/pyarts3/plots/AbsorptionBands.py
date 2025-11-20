from pyarts3.arts import AbsorptionBands, AscendingGrid, AtmPoint, AtmField, Propmat, SpeciesEnum, PropagationPathPoint
import numpy
import matplotlib
import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(data: AbsorptionBands,
         *,
         mode='normal',
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         freqs: AscendingGrid | int = 1000,
         atm: AtmPoint | AtmField | None = None,
         pol: Propmat = Propmat(1),
         species: SpeciesEnum = SpeciesEnum.AIR,
         path_point: PropagationPathPoint = PropagationPathPoint(),
         min_pm: float | None = None,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Creates a plot of the absorption by the bands in the AbsorptionBands object.

    If freqs is an index, this is the number of frequency points between the lowest
    and highest frequency in the AbsorptionBands object.  If freqs is an AscendingGrid,
    it is used directly for the x-axis.

    If atm is None, and atmospheric point is created by reading the tropical standard atmosphere
    from the ARTS database and extracting the AtmPoint at position by the path_point.  If atm is
    an AtmField, the AtmPoint at position by the path_point is extracted.

    pol is the dot product of the generated propagation matrix and defaults to the purely
    unpolarized case.

    The species is passed directly to the AbsorptionBands.spectral_propmat() method and defaults
    to all species (AIR).

    The path_point is passed directly to the AbsorptionBands.spectral_propmat() method and defaults
    to pos: [0, 0, 0], los: [0, 0].

    The mode parameter controls the style of the plot.  In 'normal' mode, all absorption lines are plotted
    as the sum of the call to AbsorptionBands.spectral_propmat().  In 'important Y X' modes, the full
    absorption is still shown, but the absorption bands object are split by X, where X is one of 'bands',
    'species', or 'isotopes'.  Y is either 'fill' or nothing.  If it is 'fill', matplotlib's
    fill_between() is used to color the regions, otherwise line plots are used to show the absorption
    of the important contributors.

    .. rubric:: Example

    .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        import numpy as np
        import pyarts3 as pyarts

        lines = pyarts.arts.AbsorptionBands.fromxml("lines/O2-66.xml")
        freq = np.linspace(20e9, 140e9, 1001)
        lines.keep_frequencies(freq[0], freq[-1])
        f, a = pyarts.plots.AbsorptionBands.plot(lines, mode="important fill bands", freqs=freq, fig=plt.figure(figsize=(16, 5)))
        a.set_yscale("log")
        a.set_xlabel("Frequency [GHz]")
        a.set_xticks(np.linspace(20, 140, 7) * 1e9, np.linspace(20, 140, 7))
        a.set_ylabel("Absorption [1/m]")
        a.legend()
        a.set_title("O$_2$ line-by-line absorption 20-140 GHz separated by bands")

    Parameters
    ----------
    data : ~pyarts3.arts.AbsorptionBands
        The AbsorptionBands object containing the data to plot.
    mode : str, optional
        The mode to use for the plot - see text for descriptions. Defaults to 'normal'.
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes] | None, optional
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
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig :
        As input if input.  Otherwise the created Figure.
    ax :
        As input if input.  Otherwise the created Axes.
    """
    if isinstance(freqs, int):
        f0 = 1e99
        f1 = -1e99
        for band in data:
            for line in data[band].lines:
                f0 = min(f0, line.f0)
                f1 = max(f1, line.f0)
        freqs = AscendingGrid(np.linspace(f0, f1, freqs))

    if atm is None:
        ws = pyarts.Workspace()
        ws.abs_speciesSet(species=[f"{band.isot}" for band in data])
        basename = "planets/Earth/afgl/tropical/"
        toa = 1 + path_point.pos[0]
        ws.atm_fieldRead(toa=toa, basename=basename, missing_is_zero=1)
        atm = ws.atm_field(*path_point.pos)
    elif isinstance(atm, AtmField):
        atm = atm(*path_point.pos)

    fig, ax = default_fig_ax(fig, ax, 1, 1, fig_kwargs={'figsize': (6, 4)})

    pm = data.spectral_propmat(
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
                'isotopes') else band.isot.spec for band in data]
            keys = list(set(keys))
            for key in keys:
                bands = data.extract_species(key)
                pm_band = bands.spectral_propmat(
                    f=freqs, atm=atm, spec=species, path_point=path_point)
                pm_bands[key] = np.einsum('ij,j->i', pm_band, pol)
        elif mode.endswith('bands'):
            keys = [band for band in data]
            for key in keys:
                bands = AbsorptionBands({key: data[key]})
                pm_band = bands.spectral_propmat(
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

        colors = matplotlib.colormaps['viridis'].resampled(len(pm_bands))

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
