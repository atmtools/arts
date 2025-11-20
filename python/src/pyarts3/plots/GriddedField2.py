""" Plotting routine for GriddedField2 """

import numpy
import matplotlib
import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.GriddedField2,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Plot a GriddedField2 as a 2D heatmap using its grids.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create a geographic field
        lats = pyarts.arts.Vector(np.linspace(-90, 90, 20))
        lons = pyarts.arts.Vector(np.linspace(-180, 175, 36))
        lon_mesh, lat_mesh = np.meshgrid(lons, lats)

        # Example: distance from equator
        data = pyarts.arts.Matrix(np.abs(lon_mesh))

        field = pyarts.arts.GriddedField2()
        field.grids = (lats, lons)
        field.data = data
        field.dataname = "Distance from Equator"

        pyarts.plots.GriddedField2.plot(field)

    Parameters
    ----------
    data : ~pyarts3.arts.GriddedField2
        A 2D gridded field with named grids
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes] | None, optional
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
    fig, ax = default_fig_ax(fig, ax, fig_kwargs={"figsize": (12, 8)})

    xgrid = data.grids[0]
    ygrid = data.grids[1]
    data = data.data
    y_mesh, x_mesh = np.meshgrid(ygrid, xgrid)

    select_flat_ax(ax, 0).pcolormesh(x_mesh, y_mesh, data, **kwargs)

    return fig, ax
