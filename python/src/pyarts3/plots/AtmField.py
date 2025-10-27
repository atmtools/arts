""" Plotting routine for profiles of the atmospheric field """

import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(
    atm_field: pyarts.arts.AtmField,
    *,
    fig=None,
    ax=None,
    alts: pyarts.arts.AscendingGrid | float = np.linspace(0, 1e5, 51),
    lats: pyarts.arts.LatGrid | float = 0,
    lons: pyarts.arts.LonGrid | float = 0,
    ygrid: pyarts.arts.Vector | None = None,
    keys: list[str] | None = None,
    **kwargs,
):
    """Plot select atmospheric field parameters by extracting a profile.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        ws = pyarts.Workspace()

        ws.atmospheric_fieldRead(toa=100e3, basename="planets/Earth/afgl/tropical/")

        pyarts.plots.AtmField.plot(ws.atmospheric_field, keys=["p", "t"])

    Parameters
    ----------
    atm_field : ~pyarts3.arts.AtmField
        An atmospheric field
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        Not used (function creates its own subplots). Accepted for API consistency.
    alts : ~pyarts3.arts.AscendingGrid | float, optional
        A grid to plot on - must after broadcast with lats and lons be 1D. Defaults to np.linspace(0, 1e5, 51).
    lats : ~pyarts3.arts.LatGrid | float, optional
        A grid to plot on - must after broadcast with alts and lons be 1D. Defaults to 0.
    lons : ~pyarts3.arts.LonGrid | float, optional
        A grid to plot on - must after broadcast with alts and lats be 1D. Defaults to 0.
    ygrid : ~pyarts3.arts.Vector | :class:`None`, optional
        Choice of y-grid for plotting.  Uses broadcasted alts if None. Defaults to None.
    keys : list, optional
        A list of keys to plot. Defaults to None for all keys in :meth:`~pyarts3.arts.AtmField.keys`.
    **kwargs : keyword arguments
        Additional keyword arguments passed to plot()

    Returns
    -------
    fig : As input
        As input.
    ax : list
        List of matplotlib axes objects.
    """
    alts, lats, lons = np.broadcast_arrays(alts, lats, lons)
    v = atm_field(alts, lats, lons)

    keys = v[0].keys() if keys is None else keys
    N = len(keys)
    n = int(np.ceil(np.sqrt(N))) + 1

    fig, ax = default_fig_ax(fig, ax, n, n, N=N, fig_kwargs={
                             'figsize': (4 * n, 4 * n), 'constrained_layout': True})

    for i in range(N):
        select_flat_ax(ax, i).plot(
            [x[keys[i]] for x in v],
            alts if ygrid is None else ygrid,
            label=keys[i],
            **kwargs
        )
    return fig, ax
