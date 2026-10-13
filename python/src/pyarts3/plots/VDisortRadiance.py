""" Plotting routine for VDisortRadiance """

import numpy
import matplotlib
import pyarts3 as pyarts
from . import DisortRadiance

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.VDisortRadiance,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         stokes: int = 0,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Plot one Stokes component of VDISORT radiance results, as :func:`pyarts3.plots.DisortRadiance.plot`.

    Pass the ``fig`` and ``ax`` of one to the other, with ``stokes=0``, to
    compare VDISORT's I with DISORT's radiance.

    Parameters
    ----------
    data : ~pyarts3.arts.VDisortRadiance
        A VDisortRadiance object containing radiance results
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes] | None, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    stokes : int, optional
        The Stokes component, 0 to 3 for I, Q, U and V. Defaults to 0.
    **kwargs : keyword arguments
        Passed to :func:`pyarts3.plots.DisortRadiance.plot` (``alt_idx``, ``azi_idx``, ``plotstyle``, ``select`` and
        the plotting options).

    Returns
    -------
    fig : matplotlib.figure.Figure
        As input if input.  Otherwise the created Figure.
    ax : matplotlib.axes.Axes
        As input if input.  Otherwise the created Axes.
    """
    component = pyarts.arts.DisortRadiance()
    component.freq_grid = data.freq_grid
    component.alt_grid = data.alt_grid
    component.azi_grid = data.azi_grid
    component.zen_grid = data.zen_grid
    component.data = pyarts.arts.Tensor4(numpy.asarray(data.data)[..., stokes])
    return DisortRadiance.plot(component, fig=fig, ax=ax, **kwargs)
