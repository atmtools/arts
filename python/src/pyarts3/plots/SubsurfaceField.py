""" Plotting routine for profiles of the subsurface field """

import numpy
import matplotlib
import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.SubsurfaceField,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         alts: pyarts.arts.AscendingGrid | float = pyarts.arts.AscendingGrid(np.linspace(-1, 0, 3)),
         lats: pyarts.arts.LatGrid | float = 0,
         lons: pyarts.arts.LonGrid | float = 0,
         ygrid: pyarts.arts.Vector | None = None,
         keys: list | None = None,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Plot select subsurface field parameters by extracting a profile.

    .. rubric:: Example

    .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        import numpy as np
        import pyarts3 as pyarts

        subsurf_field = pyarts.arts.SubsurfaceField(bottom_depth=-1)
        subsurf_field["t"] = lambda alt, lat, lon: 295 + 5 * alt * 10
        subsurf_field["rho"] = lambda alt, lat, lon: 0.977 - 0.001 * alt

        fig = plt.figure(figsize=(14, 8))
        fig, subs = pyarts.plots.SubsurfaceField.plot(subsurf_field, alts=np.linspace(-1, 0), fig=fig, keys=["t", "rho"])
        subs.flatten()[0].set_title("Temperature profile")
        subs.flatten()[1].set_title("Density profile")
        subs.flatten()[0].set_ylabel("Depth [m]")
        subs.flatten()[0].set_xlabel("Temperature [K]")
        subs.flatten()[1].set_xlabel("Density [kg/m$^3$]")

    Parameters
    ----------
    data : ~pyarts3.arts.SubsurfaceField
        A subsurface field
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes] | None, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    alts : ~pyarts3.arts.AscendingGrid | float, optional
        A grid to plot on - must after broadcast with lats and lons be 1D. Defaults to np.linspace(0, 1e5, 51).
    lats : ~pyarts3.arts.LatGrid | float, optional
        A grid to plot on - must after broadcast with alts and lons be 1D. Defaults to 0.
    lons : ~pyarts3.arts.LonGrid | float, optional
        A grid to plot on - must after broadcast with alts and lats be 1D. Defaults to 0.
    ygrid : ~pyarts3.arts.Vector | None, optional
        Choice of y-grid for plotting.  Uses broadcasted alts if None. Defaults to None.
    keys : list, optional
        A list of keys to plot. Defaults to None for all keys in :meth:`~pyarts3.arts.SubsurfaceField.keys`.
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig :
        As input if input.  Otherwise the created Figure.
    ax :
        As input if input.  Otherwise the created Axes.
    """
    alts, lats, lons = np.broadcast_arrays(alts, lats, lons)
    v = data(alts, lats, lons)

    keys = v[0].keys() if keys is None else keys
    N = len(keys)
    n = int(np.ceil(np.sqrt(N))) + 1

    fig, ax = default_fig_ax(fig, ax, n, n, N=N, fig_kwargs={'figsize': (5 * n, 5 * n)})

    for i in range(N):
        select_flat_ax(ax, i).plot(
            [x[keys[i]] for x in v],
            alts if ygrid is None else ygrid,
            label=keys[i],
            **kwargs,
        )
    return fig, ax
