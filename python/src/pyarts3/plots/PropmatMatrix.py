""" Plotting routine for PropmatMatrix """

import numpy
import matplotlib
import pyarts3 as pyarts
import numpy as np
from . import Matrix

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.PropmatMatrix,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         xgrid: pyarts.arts.Vector | None = None,
         ygrid: pyarts.arts.Vector | None = None,
         component: pyarts.arts.Propmat = pyarts.arts.Propmat([1, 0, 0, 0, 0, 0, 0]),
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """
    Plot a propagation matrix as a 2D heatmap

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        x = np.linspace(-5, 5, 50)
        y = np.linspace(-3, 3, 40)
        X, Y = np.meshgrid(x, y)
        Z = np.exp(-(X**2 + Y**2)/5)

        propmat_matrix = pyarts.arts.PropmatMatrix(np.outer(Z.flatten(), [1, 0, 0, 0, 0, 0, 0]).reshape(40, 50, 7))
        fig, ax = pyarts.plots.PropmatMatrix.plot(propmat_matrix)

    Parameters
    ----------
    data : ~pyarts3.arts.PropmatMatrix
        A matrix of propagation matrices (7-vector or 4x4 each)
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes] | None, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    xgrid : ~pyarts3.arts.Vector | None = None,
        X-axis values. If None, uses column indices. Defaults to None.
    ygrid : ~pyarts3.arts.Vector | None = None,
        Y-axis values. If None, uses row indices. Defaults to None.
    component : ~pyarts3.arts.Propmat, optional
        Choice of polarization for the heatmap.
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig :
        As input if input.  Otherwise the created Figure.
    ax :
        As input if input.  Otherwise the created Axes.
    """

    return Matrix.plot(np.einsum("ijk,k->ij", data, component),
                       fig=fig,
                       ax=ax,
                       xgrid=xgrid,
                       ygrid=ygrid,
                       **kwargs)
