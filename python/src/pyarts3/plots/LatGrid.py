""" Plotting routine for LatGrid """

import numpy
import matplotlib
import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.LatGrid,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         polar: bool = False,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
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
    data : ~pyarts3.arts.LatGrid
        A sorted grid of latitude values [-90, 90]
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes], optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    polar : bool, optional
        If True, use polar plot. Defaults to False.
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig :
        As input if input.  Otherwise the created Figure.
    ax :
        As input if input.  Otherwise the created Axes.
    """
    fig, ax = default_fig_ax(fig, ax, ax_kwargs={"subplot_kw": {'polar': polar}}, fig_kwargs={
                             'figsize': (10, 8) if polar else (10, 6)})

    if polar:
        select_flat_ax(ax, 0).plot(np.deg2rad(data), np.ones_like(data), **kwargs)
    else:
        select_flat_ax(ax, 0).plot(np.arange(len(data)), data, **kwargs)

    return fig, ax
