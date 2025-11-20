""" Plotting routine for Matrix """

import numpy
import matplotlib
import numpy as np
import pyarts3 as pyarts
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.Matrix,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         xgrid: pyarts.arts.Vector | None = None,
         ygrid: pyarts.arts.Vector | None = None,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Plot a Matrix as a 2D heatmap.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create a simple matrix
        x, y = np.meshgrid(np.linspace(-2, 2, 30), np.linspace(-2, 2, 30))
        mat = pyarts.arts.Matrix(np.exp(-(x**2 + y**2)))

        fig, ax = pyarts.plots.Matrix.plot(mat)
        ax.set_title("2D Gaussian Matrix")

    Parameters
    ----------
    data : ~pyarts3.arts.Matrix
        A 2D array of numeric values
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes], optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    xgrid : ~pyarts3.arts.Vector | None = None,
        X-axis values. If None, uses column indices. Defaults to None.
    ygrid : ~pyarts3.arts.Vector | None = None,
        Y-axis values. If None, uses row indices. Defaults to None.
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig :
        As input if input.  Otherwise the created Figure.
    ax :
        As input if input.  Otherwise the created Axes.
    """
    fig, ax = default_fig_ax(fig, ax, fig_kwargs={'figsize': (10, 8)})

    xgrid = xgrid if xgrid is not None else np.linspace(0, 1, data.shape[1])
    ygrid = ygrid if ygrid is not None else np.linspace(0, 1, data.shape[0])

    select_flat_ax(ax, 0).contourf(xgrid, ygrid, data, **kwargs)

    return fig, ax
