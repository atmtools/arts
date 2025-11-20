""" Plotting routine for SurfaceField """

import numpy
import matplotlib
import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.SurfaceField,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         lats: pyarts.arts.LatGrid = pyarts.arts.LatGrid(np.linspace(-90, 90, 50)),
         lons: pyarts.arts.LonGrid = pyarts.arts.LonGrid(np.linspace(-180, 180 * (1.0-np.finfo(float).eps), 100)),
         keys: list | None = None,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Plot surface field parameters on a grid.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        f = pyarts.arts.SurfaceField()
        f.ellipsoid = [1.0, 1]
        f['t'] = lambda lat, lon: 280 + 0.5 * lat**2 - 10 * lon
        pyarts.plots.SurfaceField.plot(f)

    Parameters
    ----------
    data : ~pyarts3.arts.SurfaceField
        A surface field containing surface parameters
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes] | None, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    lats : ~pyarts3.arts.LatGrid, optional
        Latitude grid for sampling. Defaults to [-90, 90].
    lons : ~pyarts3.arts.LonGrid, optional
        Longitude grid for sampling. Defaults to [-180, 180).
    keys : list | None, optional
        List of keys to plot. If None, plots all available keys. Defaults to None.
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig :
        As input if input.  Otherwise the created Figure.
    ax :
        As input if input.  Otherwise the created Axes.
    """
    keys = list(data.keys()) if keys is None else keys

    N = len(keys)
    n = int(np.ceil(np.sqrt(N)))

    fig, ax = default_fig_ax(fig, ax, n, n, N=N, fig_kwargs={
                             "figsize": (5 * n, 4 * n), "constrained_layout": True})

    for i, key in enumerate(keys):
        # Sample the surface field at grid points
        values = np.zeros((lats.shape[0], lons.shape[0]))
        for ii in range(lats.shape[0]):
            for jj in range(lons.shape[0]):
                point = data(lats[ii], lons[jj])
                values[ii, jj] = point[key] if key in point else np.nan

        select_flat_ax(ax, i).contourf(lons, lats, values, **kwargs)

    return fig, ax
