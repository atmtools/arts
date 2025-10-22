""" Plotting routine for SurfaceField """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    surface_field: pyarts.arts.SurfaceField,
    *,
    fig=None,
    lats: np.ndarray | None = None,
    lons: np.ndarray | None = None,
    keys: list[str] | None = None,
):
    """Plot surface field parameters on a grid.

    Parameters
    ----------
    surface_field : ~pyarts3.arts.SurfaceField
        A surface field containing surface parameters
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    lats : :class:`~numpy.ndarray` | None, optional
        Latitude grid for sampling. Defaults to None for automatic grid.
    lons : :class:`~numpy.ndarray` | None, optional
        Longitude grid for sampling. Defaults to None for automatic grid.
    keys : list[str] | None, optional
        List of keys to plot. If None, plots all available keys. Defaults to None.

    Returns
    -------
    fig : As input
        The matplotlib figure.
    axes : list
        List of matplotlib axes objects.
    """
    if keys is None:
        keys = list(surface_field.keys())
    
    if lats is None:
        lats = np.linspace(-90, 90, 50)
    if lons is None:
        lons = np.linspace(-180, 180, 100)
    
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    
    N = len(keys)
    n = int(np.ceil(np.sqrt(N)))
    
    if fig is None:
        fig = plt.figure(figsize=(5 * n, 4 * n), constrained_layout=True)
    
    axes = []
    for i, key in enumerate(keys):
        ax = fig.add_subplot(n, n, i + 1)
        
        # Sample the surface field at grid points
        values = np.zeros_like(lat_grid)
        for ii in range(lat_grid.shape[0]):
            for jj in range(lat_grid.shape[1]):
                point = surface_field(lat_grid[ii, jj], lon_grid[ii, jj])
                values[ii, jj] = point[key] if key in point else np.nan
        
        im = ax.pcolormesh(lon_grid, lat_grid, values, cmap='viridis', shading='auto')
        plt.colorbar(im, ax=ax, label=key)
        ax.set_xlabel('Longitude [°]')
        ax.set_ylabel('Latitude [°]')
        ax.set_title(key)
        axes.append(ax)
    
    return fig, axes
