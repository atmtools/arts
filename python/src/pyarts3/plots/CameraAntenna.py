""" Plotting routine for CameraAntenna """

import numpy
import matplotlib
import numpy as np
import pyarts3 as pyarts

from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.sensor.CameraAntenna,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Plot the pixel layout of a CameraAntenna in CCD coordinates.

    Each pixel is drawn at its column and row offset from the optical axis.  Rows go towards
    increasing zenith and, because the lens inverts the image, columns towards decreasing azimuth.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts

        camera = pyarts.arts.sensor.CameraAntenna.thin_lens(8, 12, 0.024, 0.036, 0.05, 1e3)

        pyarts.plots.CameraAntenna.plot(camera)

    Parameters
    ----------
    data : ~pyarts3.arts.sensor.CameraAntenna
        The camera whose pixels to plot
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes] | None, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    **kwargs : keyword arguments
        Additional keyword arguments to pass to :meth:`~matplotlib.axes.Axes.scatter`.

    Returns
    -------
    fig : matplotlib.figure.Figure
        As input if input.  Otherwise the created Figure.
    ax : matplotlib.axes.Axes
        As input if input.  Otherwise the created Axes.
    """
    fig, ax = default_fig_ax(fig, ax, fig_kwargs={'figsize': (8, 6)})
    a = select_flat_ax(ax, 0)

    rows = np.asarray(data.grid.rows, dtype=float)
    cols = np.asarray(data.grid.cols, dtype=float)
    col_grid, row_grid = np.meshgrid(cols, rows)

    a.scatter(col_grid.ravel(), row_grid.ravel(), **({"marker": "s"} | kwargs))
    a.set_xlabel("Column offset (CCD) [deg]")
    a.set_ylabel("Row offset (CCD) [deg]")
    a.set_aspect("equal")

    return fig, ax
