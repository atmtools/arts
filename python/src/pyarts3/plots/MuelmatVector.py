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
         component: None | pyarts.arts.Muelmat = None,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Plot the transmission matrix elements.

    Selects a 1D plot if the transmission matrix is 1D, otherwise a 4x4 grid of plots.

    Parameters
    ----------
    data : ~pyarts3.arts.MuelmatVector
        A vector of Mueller matrices representing the transmission properties.
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes], optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    freqs : ~pyarts3.arts.AscendingGrid, optional
        A grid of frequencies to plot. Defaults to None for no frequency grid.
    component : None | ~pyarts3.arts.Muelmat, optional
        If None, plot the dot product of each Mueller matrix with the given component.
        If provided, plot each of the 16 elements M[i,j] across the vector.
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
                             1 if component is None else 4,
                             1 if component is None else 4,
                             fig_kwargs={'figsize': (6, 4) if component is None else (12, 12)})

    if component is None:
        select_flat_ax(ax, 0).plot(freqs, np.einsum(
            "ik,k->i", data, component), **kwargs)
    else:
        for i in range(4):
            for j in range(4):
                select_flat_ax(ax, i*4 + j).plot(freqs,
                                                 data[:, i, j], **kwargs)

    return fig, ax
