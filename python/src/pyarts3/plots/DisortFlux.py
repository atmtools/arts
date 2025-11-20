""" Plotting routine for DisortFlux """

import numpy
import matplotlib
import pyarts3 as pyarts
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.DisortFlux,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         freq_idx: int | None = None,
         select: list | str = ['up', 'down_diffuse', 'down_direct'],
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Plot DISORT flux results.

    .. rubric:: Example (mock data)

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        flux = pyarts.arts.DisortFlux()
        flux.freq_grid = pyarts.arts.AscendingGrid(np.linspace(1e9, 1e12, 50))
        flux.alt_grid = pyarts.arts.DescendingGrid(np.linspace(80e3, 0, 41))
        alts_layers = 0.5 * (flux.alt_grid[:-1] + flux.alt_grid[1:])
        flux.up = pyarts.arts.Matrix(np.outer(np.linspace(10, 100, 50), np.exp(-alts_layers/20e3)))
        flux.down_diffuse = pyarts.arts.Matrix(np.outer(np.linspace(50, 150, 50), np.exp(-alts_layers/15e3)))
        flux.down_direct = pyarts.arts.Matrix(np.outer(np.linspace(200, 300, 50),  np.exp(-alts_layers/10e3)))
        fig, ax = pyarts.plots.DisortFlux.plot(flux)

    Parameters
    ----------
    data : ~pyarts3.arts.DisortFlux
        A DisortFlux object containing flux results
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes] | None, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    freq_idx : int | None, optional
        Frequency index to plot. If None, plots all frequencies. Defaults to None.
    select : list or str, optional
        Which flux components to plot: 'up', 'down_diffuse', and/or 'down_direct'.
        Defaults to all three.
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig :
        As input if input.  Otherwise the created Figure.
    ax :
        As input if input.  Otherwise the created Axes.
    """

    freq_grid = data.freq_grid
    alt_grid = data.alt_grid
    upwelling = data.up
    downwelling_diffuse = data.down_diffuse
    downwelling_direct = data.down_direct

    select = [select] if isinstance(select, str) else select
    has_up = 'up' in select
    has_down_diffuse = 'down_diffuse' in select
    has_down_direct = 'down_direct' in select
    n = sum([has_up, has_down_diffuse, has_down_direct])
    if n == 0:
        raise ValueError(
            "At least one of 'up', 'down_diffuse', or 'down_direct' must be selected.")

    fig, ax = default_fig_ax(fig, ax, 1, n, fig_kwargs={
                             'figsize': (14, 5), 'constrained_layout': True})

    idx = 0
    if freq_idx is not None:
        if has_up:
            select_flat_ax(ax, idx).plot(upwelling[freq_idx, :], alt_grid, **kwargs)
            idx += 1
        if has_down_diffuse:
            select_flat_ax(ax, idx).plot(
                downwelling_diffuse[freq_idx, :], alt_grid, **kwargs)
            idx += 1
        if has_down_direct:
            select_flat_ax(ax, idx).plot(
                downwelling_direct[freq_idx, :], alt_grid, **kwargs)
    else:
        if has_up:
            select_flat_ax(ax, idx).contourf(
                freq_grid, alt_grid[:-1], upwelling.T, **kwargs)
            idx += 1
        if has_down_diffuse:
            select_flat_ax(ax, idx).contourf(
                freq_grid, alt_grid[:-1], downwelling_diffuse.T,  **kwargs)
            idx += 1
        if has_down_direct:
            select_flat_ax(ax, idx).contourf(
                freq_grid, alt_grid[:-1], downwelling_direct.T, **kwargs)

    return fig, ax
