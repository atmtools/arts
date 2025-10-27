""" Plotting routine for LatGrid """

import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(
    grid: pyarts.arts.LatGrid,
    *,
    fig=None,
    ax=None,
    polar: bool = False,
    **kwargs
):
    """Plot a LatGrid showing latitude values.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create a latitude grid
        lats = pyarts.arts.LatGrid(np.linspace(-90, 90, 19))

        pyarts.plots.LatGrid.plot(lats, polar=True)

    Parameters
    ----------
    grid : ~pyarts3.arts.LatGrid
        A sorted grid of latitude values [-90, 90]
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    polar : bool, optional
        If True, use polar plot. Defaults to False.
    **kwargs
        Additional keyword arguments passed to plot()

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : As input
        The matplotlib axes.
    """
    fig, ax = default_fig_ax(fig, ax, ax_kwargs={"subplot_kw": {'polar': polar}}, fig_kwargs={
                             'figsize': (10, 8) if polar else (10, 6)})

    if polar:
        select_flat_ax(ax, 0).plot(np.deg2rad(grid), np.ones_like(grid), **kwargs)
    else:
        select_flat_ax(ax, 0).plot(np.arange(len(grid)), grid, **kwargs)

    return fig, ax
