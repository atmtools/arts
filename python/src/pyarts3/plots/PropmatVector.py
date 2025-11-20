""" Plotting routine for PropmatVector """

import numpy
import matplotlib
import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.PropmatVector,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         freqs: pyarts.arts.AscendingGrid | None = None,
         component: pyarts.arts.Propmat | None = None,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """
    Plot a propagation matrix.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        propmat_vec = pyarts.arts.PropmatVector(np.outer(np.exp(-np.linspace(0, 3, 100)), [1, 0, 0, 0, 0, 0, 0]))
        fig, ax = pyarts.plots.PropmatVector.plot(propmat_vec)

    Parameters
    ----------
    data : ~pyarts3.arts.PropmatVector
        A vector of propagation matrices
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes] | None, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    freqs : ~pyarts3.arts.AscendingGrid | None, optional
        Frequency or position grid for x-axis. If None, uses indices.
    component : ~pyarts3.arts.Propmat | None, optional
        If None, show grid of 16 subplots (M[i,j]). If a 7-vector, plot dot product. If int, plot compact component.
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig :
        As input if input.  Otherwise the created Figure.
    ax :
        As input if input.  Otherwise the created Axes.
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
        select_flat_ax(ax, 0).plot(freqs, np.einsum(
            "ij,j->i", data, component), **kwargs)

    return fig, ax
