""" Plotting routine for Vector """

import numpy
import matplotlib
import pyarts3 as pyarts
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.Vector,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         xgrid=None,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Plot a Vector as a line plot.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create a simple vector
        x = np.linspace(0, 2*np.pi, 50)
        vec = pyarts.arts.Vector(np.sin(x))

        pyarts.plots.Vector.plot(vec, xgrid=pyarts.arts.convert.rad2deg(x))

    Parameters
    ----------
    data : ~pyarts3.arts.Vector
        A 1D array of numeric values
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes], optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    xgrid : ~pyarts3.arts.Vector, optional
        The x-coordinates for the plot. If None, the index of the data is used.
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig :
        As input if input.  Otherwise the created Figure.
    ax :
        As input if input.  Otherwise the created Axes.
    """

    fig, ax = default_fig_ax(fig, ax, 1, 1, fig_kwargs={"figsize": (8, 6)})
    if xgrid is not None:
        select_flat_ax(ax, 0).plot(xgrid, data, **kwargs)
    else:
        select_flat_ax(ax, 0).plot(data, **kwargs)

    return fig, ax
