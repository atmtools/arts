""" Plotting routine for AziGrid """

import numpy
import matplotlib
import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.AziGrid,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         polar: bool = False,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Plot an AziGrid showing azimuth angles.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create azimuth angles (compass directions)
        azimuth = pyarts.arts.AziGrid(np.linspace(0, 360, 13)[:-1])

        fig, ax = pyarts.plots.AziGrid.plot(azimuth, polar=True)
        ax.set_xlabel("Index")
        ax.set_ylabel("Azimuth Angle [°]")
        ax.set_title("Azimuth Grid")
        ax.set_ylim(0, 360)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 1.2)
        ax.set_theta_zero_location("N")  # 0° at North (top)
        ax.set_theta_direction(-1)  # Clockwise (East = 90° clockwise from North)

    Parameters
    ----------
    data : ~pyarts3.arts.AziGrid
        A sorted grid of azimuth angles [0, 360)
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes], optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    polar : bool, optional
        If True, use polar plot. Defaults to False.
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig :
        As input if input.  Otherwise the created Figure.
    ax :
        As input if input.  Otherwise the created Axes.
    """
    fig, ax = default_fig_ax(fig, ax, ax_kwargs={"subplot_kw": {'polar': polar}}, fig_kwargs={
                             'figsize': (10, 8) if polar else (10, 6)})

    if polar:
        select_flat_ax(ax, 0).plot(np.deg2rad(data), np.ones_like(data), **kwargs)
    else:
        select_flat_ax(ax, 0).plot(np.arange(len(data)), data, **kwargs)

    return fig, ax
