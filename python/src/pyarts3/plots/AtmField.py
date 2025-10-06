""" Plotting routine for profiles of the atmospheric field """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    atm_field: pyarts.arts.AtmField,
    *,
    fig=None,
    alts: np.ndarray | float = np.linspace(0, 1e5, 51),
    lats: np.ndarray | float = 0,
    lons: np.ndarray | float = 0,
    ygrid: np.ndarray | None = None,
    keys: list[str] | None = None,
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
    alts : :class:`~numpy.ndarray` | :class:`float`, optional
        A grid to plot on - must after broadcast with lats and lons be 1D. Defaults to np.linspace(0, 1e5, 51).
    lats : :class:`~numpy.ndarray` | :class:`float`, optional
        A grid to plot on - must after broadcast with alts and lons be 1D. Defaults to 0.
    lons : :class:`~numpy.ndarray` | :class:`float`, optional
        A grid to plot on - must after broadcast with alts and lats be 1D. Defaults to 0.
    ygrid : :class:`~numpy.ndarray` | :class:`None`, optional
        Choice of y-grid for plotting.  Uses broadcasted alts if None. Defaults to None.
    keys : list, optional
        A list of keys to plot. Defaults to None for all keys in :meth:`~pyarts3.arts.AtmField.keys`.

    Returns
    -------
    fig : As input
        As input.
    subs : As input
        As input.
    """
    alts, lats, lons = np.broadcast_arrays(alts, lats, lons)
    v = atm_field(alts, lats, lons)

    keys = v[0].keys() if keys is None else keys
    N = len(keys)
    n = int(np.ceil(np.sqrt(N))) + 1

    if fig is None:
        fig = plt.figure(figsize=(4 * n, 4 * n), constrained_layout=True)

    subs = []
    for i in range(N):
        subs.append(fig.add_subplot(n, n, i + 1))
        subs[-1].plot(
            [x[keys[i]] for x in v],
            alts if ygrid is None else ygrid,
            label=keys[i],
        )
        subs[-1].legend()
    return fig, subs
