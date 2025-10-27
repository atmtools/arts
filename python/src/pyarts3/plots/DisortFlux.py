""" Plotting routine for DisortFlux """

import pyarts3 as pyarts
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(
    data: pyarts.arts.DisortFlux,
    *,
    fig=None,
    ax=None,
    freq_idx: int | None = None,
    **kwargs,
):
    """Plot DISORT flux results.

    Parameters
    ----------
    data : ~pyarts3.arts.DisortFlux
        A DisortFlux object containing flux results
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        Not used (function creates its own subplots). Accepted for API consistency.
    freq_idx : int | None, optional
        Frequency index to plot. If None, plots all frequencies. Defaults to None.
    **kwargs
        Additional keyword arguments passed to matplotlib plot/pcolormesh functions.

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : list
        List of matplotlib axes objects.
    """

    freq_grid = data.frequency_grid
    alt_grid = data.altitude_grid
    upwelling = data.up
    downwelling_diffuse = data.down_diffuse
    downwelling_direct = data.down_direct

    fig, ax = default_fig_ax(fig, ax, 1, 3, fig_kwargs={
                             'figsize': (14, 5), 'constrained_layout': True})

    if freq_idx is not None:
        select_flat_ax(ax, 0).plot(upwelling[freq_idx, :], alt_grid, **kwargs)
        select_flat_ax(ax, 1).plot(downwelling_diffuse[freq_idx, :], alt_grid, **kwargs)
        select_flat_ax(ax, 2).plot(downwelling_direct[freq_idx, :], alt_grid, **kwargs)
    else:
        extent = [freq_grid[0], freq_grid[-1], alt_grid[-1], alt_grid[0]]
        select_flat_ax(ax, 0).imshow(upwelling.T, extent=extent, **kwargs)
        select_flat_ax(ax, 1).imshow(downwelling_diffuse.T, extent=extent, **kwargs)
        select_flat_ax(ax, 2).imshow(downwelling_direct.T, extent=extent, **kwargs)

    return fig, ax
