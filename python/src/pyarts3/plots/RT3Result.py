""" Plotting routine for RT3Result """

import numpy
import matplotlib
import pyarts3 as pyarts
from .common import stream_radiance

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.rt3.RT3Result,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: numpy.ndarray | None = None,
         level: int = 0,
         azimuth: float = 0.0,
         stokes: list | None = None,
         **kwargs) -> tuple[matplotlib.figure.Figure, numpy.ndarray]:
    """Plot the RT3 radiance at one level and azimuth against the line-of-sight zenith angle.

    The radiance is summed from the solution's Fourier modes at the azimuth
    with :func:`pyarts3.arts.rt3.azimuth_radiance`.  The extra angles are
    included.  Below 90 degrees the line of sight looks up and sees the
    downward radiance, above 90 degrees it looks down and sees the upward
    radiance.  Pass the ``fig`` and ``ax`` of
    :func:`pyarts3.plots.cppvdisort.plot` or
    :func:`pyarts3.plots.RT4Result.plot` to compare the solvers.

    .. rubric:: Example

    .. code-block:: python

        import pyarts3 as pyarts

        # problem = pyarts.arts.rt3.problem_from_path(...)
        result = pyarts.arts.rt3.solve(problem)
        fig, ax = pyarts.plots.RT3Result.plot(result, level=0, azimuth=90.0, label="RT3")

    Parameters
    ----------
    data : ~pyarts3.arts.rt3.RT3Result
        The RT3 solution.
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~numpy.ndarray of ~matplotlib.axes.Axes, optional
        One Axes per Stokes component to draw on. Defaults to None for new axes.
    level : int, optional
        The level, 0 at the top. Defaults to 0.
    azimuth : float, optional
        RT3's azimuth [degree] of the propagation direction, relative to the
        beam's. Defaults to 0.
    stokes : list of int, optional
        The Stokes components to show. Defaults to all of the solution.
    **kwargs : keyword arguments
        Passed to :func:`matplotlib.axes.Axes.plot`.

    Returns
    -------
    fig : ~matplotlib.figure.Figure
        As input if input.  Otherwise the created Figure.
    ax : ~numpy.ndarray of ~matplotlib.axes.Axes
        As input if input.  Otherwise the created Axes.
    """
    phi = [numpy.radians(azimuth)]
    up = numpy.asarray(pyarts.arts.rt3.azimuth_radiance(data.up, phi))[level, 0]
    down = numpy.asarray(pyarts.arts.rt3.azimuth_radiance(data.down, phi))[level, 0]
    return stream_radiance(data.mu, up, down, fig=fig, ax=ax, stokes=stokes, **kwargs)
