""" Plotting routine for Matrix """

import pyarts3 as pyarts
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(
    matrix: pyarts.arts.Matrix,
    *,
    fig=None,
    ax=None,
    xgrid: pyarts.arts.Vector | None = None,
    ygrid: pyarts.arts.Vector | None = None,
    **kwargs
):
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
    matrix : ~pyarts3.arts.Matrix
        A 2D array of numeric values
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    xgrid : ~pyarts3.arts.Vector | None = None,
        X-axis values. If None, uses column indices. Defaults to None.
    ygrid : ~pyarts3.arts.Vector | None = None,
        Y-axis values. If None, uses row indices. Defaults to None.
    **kwargs
        Additional keyword arguments passed to imshow()

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : As input
        The matplotlib axes.
    """
    fig, ax = default_fig_ax(fig, ax, fig_kwargs={'figsize': (10, 8)})
    
    if xgrid is not None and ygrid is not None:
        select_flat_ax(ax, 0).pcolormesh(xgrid, ygrid, matrix, **kwargs)
    else:
        select_flat_ax(ax, 0).imshow(matrix, **kwargs)
    
    return fig, ax
