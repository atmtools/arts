""" Plotting routine for cppvdisort """

import numpy
import matplotlib
import pyarts3 as pyarts
from .common import stream_radiance

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.cppvdisort,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: numpy.ndarray | None = None,
         level: int = 0,
         azimuth: float = 0.0,
         stokes: list | None = None,
         **kwargs) -> tuple[matplotlib.figure.Figure, numpy.ndarray]:
    """Plot the VDISORT radiance at a layer boundary and azimuth against the line-of-sight zenith angle.

    The radiance is evaluated on the streams with ``data.u`` at the optical
    depth of the boundary.  Below 90 degrees the line of sight looks up and
    sees the downward radiance, above 90 degrees it looks down and sees the
    upward radiance.  Pass the ``fig`` and ``ax`` to
    :func:`pyarts3.plots.RT3Result.plot` or
    :func:`pyarts3.plots.RT4Result.plot` to compare the solvers.

    .. rubric:: Example

    .. code-block:: python

        import pyarts3 as pyarts

        # solution = pyarts.arts.vdisort.main_data_from_path(...)
        fig, ax = pyarts.plots.cppvdisort.plot(solution, level=0, azimuth=90.0, stokes=[0, 1], label="VDISORT")

    Parameters
    ----------
    data : ~pyarts3.arts.cppvdisort
        A solved VDISORT problem.
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~numpy.ndarray of ~matplotlib.axes.Axes, optional
        One Axes per Stokes component to draw on. Defaults to None for new axes.
    level : int, optional
        The layer boundary, 0 at the top and ``len(data.tau)`` at the bottom.
        Defaults to 0.
    azimuth : float, optional
        VDISORT's azimuth phi [degree] of the propagation direction. Defaults
        to 0.
    stokes : list of int, optional
        The Stokes components to show. Defaults to all four.
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig : ~matplotlib.figure.Figure
        As input if input.  Otherwise the created Figure.
    ax : ~numpy.ndarray of ~matplotlib.axes.Axes
        As input if input.  Otherwise the created Axes.
    """
    tau = numpy.concatenate(([0.0], numpy.asarray(data.tau)))[level]
    u = numpy.asarray(data.u(tau=[tau], phi=[numpy.radians(azimuth)]))[0, 0]
    n = len(data.weights)
    return stream_radiance(numpy.asarray(data.mu)[:n], u[:n], u[n:], fig=fig, ax=ax, stokes=stokes, **kwargs)
