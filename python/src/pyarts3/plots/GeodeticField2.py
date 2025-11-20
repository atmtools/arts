""" Plotting routine for GeodeticField2 """

import numpy
import matplotlib
import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.GeodeticField2,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Plot a GeodeticField2 as a map on latitude/longitude grid.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create a geographic field
        lats = pyarts.arts.LatGrid(np.linspace(-90, 90, 20))
        lons = pyarts.arts.LonGrid(np.linspace(-180, 175, 36))
        lon_mesh, lat_mesh = np.meshgrid(lons, lats)

        # Example: distance from equator
        data = pyarts.arts.Matrix(np.abs(lat_mesh))

        field = pyarts.arts.GeodeticField2()
        field.grids = (lats, lons)
        field.data = data
        field.dataname = "Distance from Equator"

        pyarts.plots.GeodeticField2.plot(field)

    Parameters
    ----------
    data : ~pyarts3.arts.GeodeticField2
        A 2D geodetic field with lat/lon grids
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
    fig, ax = default_fig_ax(fig, ax, fig_kwargs={"figsize": (12, 8)})

    lat_grid = data.grids[0]
    lon_grid = data.grids[1]
    data = data.data
    lon_mesh, lat_mesh = np.meshgrid(lon_grid, lat_grid)

    select_flat_ax(ax, 0).pcolormesh(lon_mesh, lat_mesh, data, **kwargs)

    return fig, ax
