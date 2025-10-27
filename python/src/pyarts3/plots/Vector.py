""" Plotting routine for Vector """

import pyarts3 as pyarts
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(
    vector: pyarts.arts.Vector,
    *,
    fig=None,
    ax=None,
    **kwargs
):
    """Plot a Vector as a line plot.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create a simple vector
        vec = pyarts.arts.Vector(np.sin(np.linspace(0, 2*np.pi, 50)))

        pyarts.plots.Vector.plot(vec)

    Parameters
    ----------
    vector : ~pyarts3.arts.Vector
        A 1D array of numeric values
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    **kwargs
        Additional keyword arguments passed to plot()

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : As input
        The matplotlib axes.
    """

    fig, ax = default_fig_ax(fig, ax, 1, 1, fig_kwargs={"figsize": (8, 6)})

    select_flat_ax(ax, 0).plot(vector, **kwargs)

    return fig, ax
