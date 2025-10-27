""" Plotting routine for LonGrid """

import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(
    data: pyarts.arts.LonGrid,
    *,
    fig=None,
    ax=None,
    polar: bool = False,
    **kwargs
):
    """Plot a LonGrid showing longitude values.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create a longitude grid
        lons = pyarts.arts.LonGrid(np.linspace(-180, 175, 36))

        pyarts.plots.LonGrid.plot(lons, polar=True)

    Parameters
    ----------
    data : ~pyarts3.arts.LonGrid
        A sorted grid of longitude values [-180, 180)
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
        select_flat_ax(ax, 0).plot(np.deg2rad(data), np.ones_like(data), **kwargs)
    else:
        select_flat_ax(ax, 0).plot(np.arange(len(data)), data, **kwargs)

    return fig, ax
