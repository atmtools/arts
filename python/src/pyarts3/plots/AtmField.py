""" Plotting routine for profiles of the atmospheric field """

import numpy
import matplotlib
import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def vmr_scale(unit):
    scale = {
        'm': '‰',
        'µ': 'ppmv',
        'n': 'ppbv',
        'p': 'pptv',
        'f': 'ppqv'
    }

    return scale.get(unit, f"{unit}VMR")


def get_label(x, unit=''):
    if isinstance(x, pyarts.arts.AtmKey):
        match x:
            case pyarts.arts.AtmKey.t: return f"Temperature [{unit}K]"
            case pyarts.arts.AtmKey.p: return f"Pressure [{unit}Pa]"
            case pyarts.arts.AtmKey.wind_u: return f"Wind field U [{unit}m/s]"
            case pyarts.arts.AtmKey.wind_v: return f"Wind field V [{unit}m/s]"
            case pyarts.arts.AtmKey.wind_w: return f"Wind field W [{unit}m/s]"
            case pyarts.arts.AtmKey.mag_u: return f"Magnetic field U [{unit}T]"
            case pyarts.arts.AtmKey.mag_v: return f"Magnetic field V [{unit}T]"
            case pyarts.arts.AtmKey.mag_w: return f"Magnetic field W [{unit}T]"

    if isinstance(x, pyarts.arts.SpeciesEnum) or isinstance(x, pyarts.arts.SpeciesIsotope):
        return f"{x} [{vmr_scale(unit)}]"

    return f"{x} [{unit if unit != '' else '-'}]"


def natural_scale(x):
    x = np.array(x)
    ix = np.nonzero(x)[0]
    if len(ix) == 0:
        return '', x

    xmin = x[ix].min()
    [scl0, xmin0] = pyarts.arts.convert.metric_prefix(xmin)
    scl0 = scl0 if scl0 != 'u' else 'µ'
    scl0 = '' if scl0 == '' else scl0

    xmax = x[ix].max()
    [scl1, xmax0] = pyarts.arts.convert.metric_prefix(xmax)
    scl1 = scl1 if scl1 != 'u' else 'µ'
    scl1 = '' if scl1 == '' else scl1

    if (scl0 == scl1) or (xmin < xmin0 and xmax < xmax0) or (xmin > xmin0 and xmax > xmax0):
        return scl1, x / xmax * xmax0

    return '', x


def plot(data: pyarts.arts.AtmField,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         alts: pyarts.arts.AscendingGrid | float = np.linspace(0, 1e5, 51),
         lats: pyarts.arts.LatGrid | float = 0,
         lons: pyarts.arts.LonGrid | float = 0,
         ygrid: pyarts.arts.Vector | None = None,
         keys: list[str] | None = None,
         apply_natural_scale: bool = False,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Plot select atmospheric field parameters by extracting a profile.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        ws = pyarts.Workspace()

        ws.atm_fieldRead(toa=100e3, basename="planets/Earth/afgl/tropical/")

        pyarts.plots.AtmField.plot(ws.atm_field, keys=["p", "t"])

    Parameters
    ----------
    data : ~pyarts3.arts.AtmField
        An atmospheric field
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes], optional
        The matplotlib axes to draw on. Defaults to None for new axes.
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
    apply_natural_scale : bool, optional
        Whether to apply natural scaling to each parameter. Defaults to False.
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

    fig, ax = default_fig_ax(fig, ax, n, n, N=N, fig_kwargs={
                             'figsize': (4 * n, 4 * n), 'constrained_layout': True})

    for i in range(N):
        if apply_natural_scale:
            unit, x = natural_scale([x[keys[i]] for x in v])
        else:
            unit = ''
            x = [x[keys[i]] for x in v]
        select_flat_ax(ax, i).plot(
            x,
            alts if ygrid is None else ygrid,
            label=get_label(keys[i], unit),
            **kwargs
        )
    return fig, ax
