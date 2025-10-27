""" Plotting routine for PropmatMatrix """

import pyarts3 as pyarts
import numpy as np
from . import Matrix

__all__ = [
    'plot',
]


def plot(
    data: pyarts.arts.PropmatMatrix,
    *,
    fig=None,
    ax=None,
    xgrid: pyarts.arts.Vector | None = None,
    ygrid: pyarts.arts.Vector | None = None,
    component: pyarts.arts.Propmat = pyarts.arts.Propmat([1, 0, 0, 0, 0, 0, 0]),
    **kwargs
):
    """
    Plot a propagation matrix as a 2D heatmap

    Parameters
    ----------
    data : ~pyarts3.arts.PropmatMatrix
        A matrix of propagation matrices (7-vector or 4x4 each)
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes or list, optional
        Axes to plot on. If None, new axes are created.
    xgrid : ~pyarts3.arts.Vector | None = None,
        X-axis values. If None, uses column indices. Defaults to None.
    ygrid : ~pyarts3.arts.Vector | None = None,
        Y-axis values. If None, uses row indices. Defaults to None.
    component : ~pyarts3.arts.Propmat, optional
        Choice of polarization for the heatmap.
    **kwargs
        Additional keyword arguments passed to imshow().

    Returns
    -------
    fig : Figure
        The matplotlib figure.
    ax : list or Axes
        List of axes (grid mode) or single axes (dot product/element/component mode).
    """

    return Matrix.plot(np.einsum("ijk,k->ij", data, component),
                       fig=fig,
                       ax=ax,
                       xgrid=xgrid,
                       ygrid=ygrid,
                       **kwargs)
