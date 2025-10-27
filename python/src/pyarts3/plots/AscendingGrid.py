""" Plotting routine for AscendingGrid """

import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(
    grid: pyarts.arts.AscendingGrid,
    *,
    fig=None,
    ax=None,
    **kwargs
):
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
    grid : ~pyarts3.arts.AscendingGrid
        A sorted ascending grid of values
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    **kwargs : keyword arguments
        Additional keyword arguments passed to plot()

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : As input
        The matplotlib axes.
    """
    fig, ax = default_fig_ax(fig, ax, 1, 1, fig_kwargs={
                             'figsize': (10, 6), 'constrained_layout': True})

    indices = np.arange(len(grid))
    select_flat_ax(ax, 0).plot(indices, grid, **kwargs)

    return fig, ax
