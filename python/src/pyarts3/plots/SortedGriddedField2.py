""" Plotting routine for SortedGriddedField2 """

import numpy
import matplotlib
import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.SortedGriddedField2,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Plot a SortedGriddedField2 as a 2D heatmap using its sorted grids.

    Parameters
    ----------
    data : ~pyarts3.arts.SortedGriddedField2
        A 2D sorted gridded field with ascending grids
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes], optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig :
        As input if input.  Otherwise the created Figure.
    ax :
        As input if input.  Otherwise the created Axes.
    """
    fig, ax = default_fig_ax(fig, ax, fig_kwargs={'figsize': (10, 6)})

    # Get grids and data
    xgrid = data.grids[0]
    ygrid = data.grids[1]
    data = data.data
    y_mesh, x_mesh = np.meshgrid(ygrid, xgrid)

    select_flat_ax(ax, 0).pcolormesh(x_mesh, y_mesh, data, **kwargs)

    return fig, ax
