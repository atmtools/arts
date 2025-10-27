""" Plotting routine for StokvecVector """

import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(
    data: pyarts.arts.StokvecVector,
    *,
    fig=None,
    ax=None,
    freqs: np.ndarray | None = None,
    component: pyarts.arts.Stokvec | None = None,
    **kwargs
):
    """
    Plot the Stokes vectors.

    Parameters
    ----------
    data : ~pyarts3.arts.StokvecVector
        A vector of Stokes vectors (each with 4 components: I, Q, U, V)
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes or list, optional
        Axes to plot on. If None, new axes are created.
    freqs : :class:`~numpy.ndarray` | None, optional
        Frequency or position grid for x-axis. If None, uses indices.
    component : ~pyarts3.arts.Stokvec | None, optional
        If None, show grid of 4 subplots (I,Q,U,V). If a 4-vector, plot dot product with each sample.
    **kwargs : keyword arguments
        Additional keyword arguments passed to matplotlib plot().

    Returns
    -------
    fig : Figure
        The matplotlib figure.
    ax : list or Axes
        List of axes (grid mode) or single axes (dot product mode).
    """

    freqs = np.arange(data.shape[0]) if freqs is None else freqs

    if component is None:
        fig, ax = default_fig_ax(fig, ax, 4, 1, fig_kwargs={
                                 'figsize': (8, 24), 'constrained_layout': True})
        select_flat_ax(ax, 0).plot(freqs, data[:, 0], label="I", **kwargs)
        select_flat_ax(ax, 1).plot(freqs, data[:, 1], label="Q", **kwargs)
        select_flat_ax(ax, 2).plot(freqs, data[:, 2], label="U", **kwargs)
        select_flat_ax(ax, 3).plot(freqs, data[:, 3], label="V", **kwargs)
    else:
        component = pyarts.arts.Stokvec(component)
        fig, ax = default_fig_ax(fig, ax, 1, 1, fig_kwargs={'figsize': (8, 6)})
        select_flat_ax(ax, 0).plot(freqs, np.einsum(
            'ij,j->i', data, component), **kwargs)

    return fig, ax
