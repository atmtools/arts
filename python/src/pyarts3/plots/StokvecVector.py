""" Plotting routine for StokvecVector """

import numpy
import matplotlib
import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.StokvecVector,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         freqs: pyarts.arts.AscendingGrid | None = None,
         component: pyarts.arts.Stokvec | None = None,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """
    Plot the Stokes vectors.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        stokvec_vec = pyarts.arts.StokvecVector(np.outer(np.sin(np.linspace(0, 3, 100)), [1, 0.5, 0.3, 0.1]))
        fig, ax = pyarts.plots.StokvecVector.plot(stokvec_vec)

    Parameters
    ----------
    data : ~pyarts3.arts.StokvecVector
        A vector of Stokes vectors (each with 4 components: I, Q, U, V)
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes] | None, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    freqs : ~pyarts3.arts.AscendingGrid | None, optional
        Frequency or position grid for x-axis. If None, uses indices.
    component : ~pyarts3.arts.Stokvec | None, optional
        If None, show grid of 4 subplots (I,Q,U,V). If a 4-vector, plot dot product with each sample.
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig :
        As input if input.  Otherwise the created Figure.
    ax :
        As input if input.  Otherwise the created Axes.
    """

    freqs = np.arange(data.shape[0]) if freqs is None else freqs

    if component is None:
        fig, ax = default_fig_ax(fig, ax, 2, 2, fig_kwargs={
                                 'figsize': (12, 12), 'constrained_layout': True})
        select_flat_ax(ax, 0).plot(freqs, data[:, 0], **kwargs)
        select_flat_ax(ax, 1).plot(freqs, data[:, 1], **kwargs)
        select_flat_ax(ax, 2).plot(freqs, data[:, 2], **kwargs)
        select_flat_ax(ax, 3).plot(freqs, data[:, 3], **kwargs)
    else:
        component = pyarts.arts.Stokvec(component)
        fig, ax = default_fig_ax(fig, ax, 1, 1, fig_kwargs={'figsize': (8, 6)})
        select_flat_ax(ax, 0).plot(freqs, np.einsum(
            'ij,j->i', data, component), **kwargs)

    return fig, ax
