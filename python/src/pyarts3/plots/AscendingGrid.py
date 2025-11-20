""" Plotting routine for AscendingGrid """

import numpy
import matplotlib
import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.AscendingGrid,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Plot an AscendingGrid as a scatter/line plot.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create a frequency grid
        freqs = pyarts.arts.AscendingGrid(np.logspace(9, 12, 20))

        pyarts.plots.AscendingGrid.plot(freqs)

    Parameters
    ----------
    data : ~pyarts3.arts.AscendingGrid
        A sorted ascending grid of values
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes], optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig :
        As input if input.  Otherwise the created Figure.
    ax :
        As input if input.  Otherwise the created Axes.
    """
    fig, ax = default_fig_ax(fig, ax, 1, 1, fig_kwargs={
                             'figsize': (10, 6), 'constrained_layout': True})

    indices = np.arange(len(data))
    select_flat_ax(ax, 0).plot(indices, data, **kwargs)

    return fig, ax
