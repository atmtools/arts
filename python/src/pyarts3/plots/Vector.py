""" Plotting routine for Vector """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    vector: pyarts.arts.Vector,
    *,
    fig=None,
    ax=None,
    xgrid: np.ndarray | None = None,
    xlabel: str = "Index",
    ylabel: str = "Value",
    title: str = "Vector",
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

        pyarts.plots.Vector.plot(vec, xlabel="Sample", ylabel="Amplitude")

    Parameters
    ----------
    vector : ~pyarts3.arts.Vector
        A 1D array of numeric values
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    xgrid : :class:`~numpy.ndarray` | None, optional
        X-axis values. If None, uses indices. Defaults to None.
    xlabel : str, optional
        Label for x-axis. Defaults to "Index".
    ylabel : str, optional
        Label for y-axis. Defaults to "Value".
    title : str, optional
        Plot title. Defaults to "Vector".
    **kwargs
        Additional keyword arguments passed to plot()

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : As input
        The matplotlib axes.
    """
    if fig is None:
        fig = plt.figure(figsize=(8, 6))
    
    if ax is None:
        ax = fig.add_subplot(1, 1, 1)
    
    if xgrid is None:
        xgrid = np.arange(len(vector))
    
    ax.plot(xgrid, vector, **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    
    return fig, ax
