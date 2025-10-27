""" Plotting routine for transmission matrix """

import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(
    transmission_matrix: pyarts.arts.MuelmatVector,
    *,
    fig=None,
    ax=None,
    freqs: pyarts.arts.AscendingGrid = None,
    component: None | pyarts.arts.Muelmat = None,
    **kwargs,
):
    """Plot the transmission matrix elements.

    Selects a 1D plot if the transmission matrix is 1D, otherwise a 4x4 grid of plots.

    Parameters
    ----------
    transmission_matrix : ~pyarts3.arts.MuelmatVector
        A vector of Mueller matrices representing the transmission properties.
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes or list, optional
        The matplotlib axes to draw on. For 1D plot: single axes. For 4x4 plot: 16 axes.
        Defaults to None for new axes.
    freqs : ~pyarts3.arts.AscendingGrid, optional
        A grid of frequencies to plot. Defaults to None for no frequency grid.
    component : None | ~pyarts3.arts.Muelmat, optional
        If None, plot the dot product of each Mueller matrix with the given component.
        If provided, plot each of the 16 elements M[i,j] across the vector.
    **kwargs : optional
        Additional keyword arguments passed to matplotlib plot().

    Returns
    -------
    fig : As input
        As input.
    ax : As input
        As input. Single axis for 1D plot, or list of 16 axes for 4x4 plot.
    """
    freqs = np.arange(len(transmission_matrix)) if freqs is None else freqs

    fig, ax = default_fig_ax(fig, ax,
                             1 if component is None else 4,
                             1 if component is None else 4,
                             fig_kwargs={'figsize': (6, 4) if component is None else (12, 12)})

    if component is None:
        select_flat_ax(ax, 0).plot(freqs, np.einsum(
            "ijk,jk->i", transmission_matrix, component), **kwargs)
    else:
        for i in range(4):
            for j in range(4):
                select_flat_ax(ax, i*4 + j).plot(freqs,
                                                 transmission_matrix[:, i, j], **kwargs)

    return fig, ax
