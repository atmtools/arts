""" Plotting routine for profiles of the subsurface field """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    subsurf_field: pyarts.arts.SubsurfaceField,
    *,
    fig=None,
    ax=None,
    alts: np.ndarray | float = np.linspace(-1, 0, 3),
    lats: np.ndarray | float = 0,
    lons: np.ndarray | float = 0,
    ygrid: np.ndarray | None = None,
    keys: list[str] | None = None,
    **kwargs,
):
    """Plot select subsurface field parameters by extracting a profile.

    Parameters
    ----------
    subsurf_field : ~pyarts3.arts.SubsurfaceField
        A subsurface field
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes or list of Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    alts : :class:`~numpy.ndarray` | :class:`float`, optional
        A grid to plot on - must after broadcast with lats and lons be 1D. Defaults to np.linspace(0, 1e5, 51).
    lats : :class:`~numpy.ndarray` | :class:`float`, optional
        A grid to plot on - must after broadcast with alts and lons be 1D. Defaults to 0.
    lons : :class:`~numpy.ndarray` | :class:`float`, optional
        A grid to plot on - must after broadcast with alts and lats be 1D. Defaults to 0.
    ygrid : :class:`~numpy.ndarray` | :class:`None`, optional
        Choice of y-grid for plotting.  Uses broadcasted alts if None. Defaults to None.
    keys : list, optional
        A list of keys to plot. Defaults to None for all keys in :meth:`~pyarts3.arts.SubsurfaceField.keys`.

    Returns
    -------
    fig : As input
        As input.
    ax : list
        List of matplotlib axes objects.
    """
    alts, lats, lons = np.broadcast_arrays(alts, lats, lons)
    v = subsurf_field(alts, lats, lons)

    keys = v[0].keys() if keys is None else keys
    N = len(keys)
    n = int(np.ceil(np.sqrt(N))) + 1

    if fig is None:
        fig = plt.figure(figsize=(5 * n, 5 * n))

    if ax is None:
        ax = []
        for i in range(N):
            ax.append(fig.add_subplot(n, n, i + 1))
    
    for i in range(N):
        ax[i].plot(
            [x[keys[i]] for x in v],
            alts if ygrid is None else ygrid,
            label=keys[i],
            **kwargs,
        )
        ax[i].legend()
    return fig, ax
