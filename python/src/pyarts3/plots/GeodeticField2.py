""" Plotting routine for GeodeticField2 """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    gridded_field: pyarts.arts.GeodeticField2,
    *,
    fig=None,
    ax=None,
    xlabel: str = "Longitude [°]",
    ylabel: str = "Latitude [°]",
    title: str | None = None,
    colorbar: bool = True,
    cmap: str = "viridis",
    **kwargs
):
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
    gridded_field : ~pyarts3.arts.GeodeticField2
        A 2D geodetic field with lat/lon grids
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    xlabel : str, optional
        Label for x-axis. Defaults to "Longitude [°]".
    ylabel : str, optional
        Label for y-axis. Defaults to "Latitude [°]".
    title : str | None, optional
        Plot title. If None, uses data name. Defaults to None.
    colorbar : bool, optional
        Whether to show colorbar. Defaults to True.
    cmap : str, optional
        Colormap name. Defaults to "viridis".
    **kwargs
        Additional keyword arguments passed to pcolormesh()

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : As input
        The matplotlib axes.
    """
    if fig is None:
        fig = plt.figure(figsize=(12, 8))
    
    if ax is None:
        ax = fig.add_subplot(1, 1, 1)
    
    # Get grids and data - assuming grids[0] is lat, grids[1] is lon
    lat_grid = gridded_field.grids[0]
    lon_grid = gridded_field.grids[1]
    data = gridded_field.data
    
    # Create mesh
    lon_mesh, lat_mesh = np.meshgrid(lon_grid, lat_grid)
    
    # Plot using pcolormesh
    im = ax.pcolormesh(lon_mesh, lat_mesh, data, cmap=cmap, shading='auto', **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    if colorbar:
        plt.colorbar(im, ax=ax, label=gridded_field.dataname if title is None else title)
    
    ax.set_title(title if title is not None else gridded_field.dataname)
    
    return fig, ax
