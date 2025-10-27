""" Plotting routine for SurfaceField """

import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(
    data: pyarts.arts.SurfaceField,
    *,
    fig=None,
    ax=None,
    lats: pyarts.arts.LatGrid | None = None,
    lons: pyarts.arts.LonGrid | None = None,
    keys: list[str] | None = None,
    **kwargs,
):
    """Plot surface field parameters on a grid.

    Parameters
    ----------
    data : ~pyarts3.arts.SurfaceField
        A surface field containing surface parameters
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        Not used (function creates its own subplots). Accepted for API consistency.
    lats : ~pyarts3.arts.LatGrid | None, optional
        Latitude grid for sampling. Defaults to None for automatic grid.
    lons : ~pyarts3.arts.LonGrid | None, optional
        Longitude grid for sampling. Defaults to None for automatic grid.
    keys : list[str] | None, optional
        List of keys to plot. If None, plots all available keys. Defaults to None.
    **kwargs
        Additional keyword arguments passed to matplotlib pcolormesh function.

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : list
        List of matplotlib axes objects.
    """
    keys = list(data.keys()) if keys is None else keys
    lats = pyarts.arts.LatGrid(np.linspace(-90, 90, 50)) if lats is None else lats
    lons = pyarts.arts.LonGrid(np.linspace(-180, 180, 100)) if lons is None else lons
    lon_grid, lat_grid = np.meshgrid(lons, lats)

    N = len(keys)
    n = int(np.ceil(np.sqrt(N)))

    fig, ax = default_fig_ax(fig, ax, n, n, N=N, fig_kwargs={
                             "figsize": (5 * n, 4 * n), "constrained_layout": True})

    for i, key in enumerate(keys):
        # Sample the surface field at grid points
        values = np.zeros_like(lat_grid)
        for ii in range(lat_grid.shape[0]):
            for jj in range(lat_grid.shape[1]):
                point = data(lat_grid[ii, jj], lon_grid[ii, jj])
                values[ii, jj] = point[key] if key in point else np.nan

        select_flat_ax(ax, i).pcolormesh(lon_grid, lat_grid, values, **kwargs)

    return fig, ax
