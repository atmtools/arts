""" Plotting routine for transmission matrix """

import numpy
import matplotlib
import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.MuelmatVector,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         freqs: pyarts.arts.AscendingGrid = None,
         component: pyarts.arts.Muelmat | None  = None,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Plot the transmission matrix elements.

    Selects a 1D plot if the transmission matrix is 1D, otherwise a 4x4 grid of plots.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        muelmat_vec = pyarts.arts.MuelmatVector(np.cos(np.linspace(0, 3, 1000)))
        fig, ax = pyarts.plots.MuelmatVector.plot(muelmat_vec)

    Parameters
    ----------
    data : ~pyarts3.arts.MuelmatVector
        A vector of Mueller matrices representing the transmission properties.
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes] | None, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    freqs : ~pyarts3.arts.AscendingGrid, optional
        A grid of frequencies to plot. Defaults to None for no frequency grid.
    component : pyarts.arts.Muelmat | None, optional
        If None, plot all 16 elements M[i,j] of the Mueller matrix across the vector.
        If provided, plot the dot product of each Mueller matrix with the given component.
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig :
        As input if input.  Otherwise the created Figure.
    ax :
        As input if input.  Otherwise the created Axes.
    """
    freqs = np.arange(len(data)) if freqs is None else freqs

    fig, ax = default_fig_ax(fig, ax,
                             4 if component is None else 1,
                             4 if component is None else 1,
                             fig_kwargs={'figsize': (12, 12) if component is None else (6, 4)})

    if component is None:
        for i in range(4):
            for j in range(4):
                a = select_flat_ax(ax, i*4 + j)
                a.plot(freqs, data[:, i, j], **kwargs)
    else:
        r = np.einsum("ijk,jk->i", data, component)
        select_flat_ax(ax, 0).plot(freqs, r, **kwargs)

    return fig, ax
