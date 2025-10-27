""" Plotting routine for PropmatVector """

import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(
    data: pyarts.arts.PropmatVector,
    *,
    fig=None,
    ax=None,
    freqs: pyarts.arts.AscendingGrid | None = None,
    component: pyarts.arts.Propmat | None = None,
    **kwargs
):
    """
    Plot a propagation matrix.

    Parameters
    ----------
    data : ~pyarts3.arts.PropmatVector
        A vector of propagation matrices
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes or list, optional
        Axes to plot on. If None, new axes are created.
    freqs : ~pyarts3.arts.AscendingGrid or None, optional
        Frequency or position grid for x-axis. If None, uses indices.
    component : ~pyarts3.arts.Propmat or None, optional
        If None, show grid of 16 subplots (M[i,j]). If a 7-vector, plot dot product. If int, plot compact component.
    **kwargs
        Additional keyword arguments passed to matplotlib plot().

    Returns
    -------
    fig : Figure
        The matplotlib figure.
    ax : list or Axes
        List of axes (grid mode) or single axes (dot product/element/component mode).
    """
    fig, ax = default_fig_ax(fig, ax,
                             1 if component is not None else 4,
                             1 if component is not None else 4,
                             fig_kwargs={'figsize': (10, 6) if component is not None else (12, 12), 'constrained_layout': True})

    freqs = np.arange(len(data)) if freqs is None else freqs

    if component is None:
        A = data[:, 0]
        B = data[:, 1]
        C = data[:, 2]
        D = data[:, 3]
        U = data[:, 4]
        V = data[:, 5]
        W = data[:, 6]
        M = [A,  B,  C,  D,
             B,  A,  U,  V,
             C, -U,  A,  W,
             D, -V, -W,  A]
        for i in range(16):
            select_flat_ax(ax, i).plot(freqs, M[i], **kwargs)
    else:
        print (data.shape, component.shape)
        select_flat_ax(ax, 0).plot(freqs, np.einsum(
            "ij,j->i", data, component), **kwargs)

    return fig, ax
