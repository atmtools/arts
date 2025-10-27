""" Plotting routine for GriddedField2 """

import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(
    data: pyarts.arts.GriddedField2,
    *,
    fig=None,
    ax=None,
    **kwargs
):
    """Plot a GriddedField2 as a 2D heatmap using its grids.

    Parameters
    ----------
    data : ~pyarts3.arts.GriddedField2
        A 2D gridded field with named grids
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    **kwargs
        Additional keyword arguments passed to pcolormesh()

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : As input
        The matplotlib axes.
    """
    fig, ax = default_fig_ax(fig, ax, fig_kwargs={"figsize": (12, 8)})

    xgrid = data.grids[0]
    ygrid = data.grids[1]
    data = data.data
    y_mesh, x_mesh = np.meshgrid(ygrid, xgrid)

    select_flat_ax(ax, 0).pcolormesh(x_mesh, y_mesh, data, **kwargs)

    return fig, ax
