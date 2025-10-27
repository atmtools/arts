""" Plotting routine for AzimuthGrid """

import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(
    grid: pyarts.arts.AzimuthGrid,
    *,
    fig=None,
    ax=None,
    polar: bool = False,
    **kwargs
):
    """Plot an AzimuthGrid showing azimuth angles.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create azimuth angles (compass directions)
        azimuth = pyarts.arts.AzimuthGrid(np.linspace(0, 360, 13)[:-1])

        fig, ax = pyarts.plots.AzimuthGrid.plot(azimuth, polar=True)
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
    grid : ~pyarts3.arts.AzimuthGrid
        A sorted grid of azimuth angles [0, 360)
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    polar : bool, optional
        If True, use polar plot. Defaults to False.
    **kwargs
        Additional keyword arguments passed to plot()

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : As input
        The matplotlib axes.
    """
    fig, ax = default_fig_ax(fig, ax, ax_kwargs={"subplot_kw": {'polar': polar}}, fig_kwargs={
                             'figsize': (10, 8) if polar else (10, 6)})

    if polar:
        select_flat_ax(ax, 0).plot(np.deg2rad(grid), np.ones_like(grid), **kwargs)
    else:
        select_flat_ax(ax, 0).plot(np.arange(len(grid)), grid, **kwargs)

    return fig, ax
