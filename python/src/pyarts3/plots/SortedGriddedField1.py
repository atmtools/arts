""" Plotting routine for SortedGriddedField1 """

import numpy
import matplotlib
import pyarts3 as pyarts
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.SortedGriddedField1,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Plot a SortedGriddedField1 as a line plot using its grid.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create a gridded field
        freqs = pyarts.arts.AscendingGrid(np.linspace(1e9, 1e12, 100))
        data = pyarts.arts.Vector(np.exp(-((freqs - 5e11)/1e11)**2))

        field = pyarts.arts.SortedGriddedField1()
        field.grids = (freqs,)
        field.data = data
        field.gridnames = ("Frequency [Hz]",)
        field.dataname = "Response"

        pyarts.plots.SortedGriddedField1.plot(field)

    Parameters
    ----------
    data : ~pyarts3.arts.SortedGriddedField1
        A 1D sorted gridded field with named grid
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
    fig, ax = default_fig_ax(fig, ax, fig_kwargs={'figsize': (10, 6)})

    # Get grid and data
    xgrid = data.grids[0]
    y = data.data

    select_flat_ax(ax, 0).plot(xgrid, y, label=data.dataname, **kwargs)

    return fig, ax
