""" Plotting routine for Matrix """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    matrix: pyarts.arts.Matrix,
    *,
    fig=None,
    ax=None,
    xgrid: np.ndarray | None = None,
    ygrid: np.ndarray | None = None,
    xlabel: str = "Column",
    ylabel: str = "Row",
    title: str = "Matrix",
    colorbar: bool = True,
    cmap: str = "viridis",
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

        pyarts.plots.Matrix.plot(mat, title="2D Gaussian")

    Parameters
    ----------
    matrix : ~pyarts3.arts.Matrix
        A 2D array of numeric values
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    xgrid : :class:`~numpy.ndarray` | None, optional
        X-axis values. If None, uses column indices. Defaults to None.
    ygrid : :class:`~numpy.ndarray` | None, optional
        Y-axis values. If None, uses row indices. Defaults to None.
    xlabel : str, optional
        Label for x-axis. Defaults to "Column".
    ylabel : str, optional
        Label for y-axis. Defaults to "Row".
    title : str, optional
        Plot title. Defaults to "Matrix".
    colorbar : bool, optional
        Whether to show colorbar. Defaults to True.
    cmap : str, optional
        Colormap name. Defaults to "viridis".
    **kwargs
        Additional keyword arguments passed to imshow()

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : As input
        The matplotlib axes.
    """
    if fig is None:
        fig = plt.figure(figsize=(10, 8))
    
    if ax is None:
        ax = fig.add_subplot(1, 1, 1)
    
    # Use pcolormesh if grids are provided, otherwise use imshow
    if xgrid is not None and ygrid is not None:
        im = ax.pcolormesh(xgrid, ygrid, matrix, cmap=cmap, **kwargs)
    else:
        im = ax.imshow(matrix, aspect='auto', cmap=cmap, origin='lower', **kwargs)
    
    if colorbar:
        plt.colorbar(im, ax=ax)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    
    return fig, ax
