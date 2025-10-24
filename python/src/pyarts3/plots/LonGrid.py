""" Plotting routine for LonGrid """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    grid: pyarts.arts.LonGrid,
    *,
    fig=None,
    ax=None,
    polar: bool = False,
    xlabel: str = "Index",
    ylabel: str = "Longitude [°]",
    title: str = "Longitude Grid",
    marker: str = 'o',
    **kwargs
):
    """Plot a LonGrid showing longitude values.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create a longitude grid
        lons = pyarts.arts.LonGrid(np.linspace(-180, 175, 36))

        pyarts.plots.LonGrid.plot(lons, polar=True)

    Parameters
    ----------
    grid : ~pyarts3.arts.LonGrid
        A sorted grid of longitude values [-180, 180)
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    polar : bool, optional
        If True, use polar plot. Defaults to False.
    xlabel : str, optional
        Label for x-axis (not used in polar). Defaults to "Index".
    ylabel : str, optional
        Label for y-axis (not used in polar). Defaults to "Longitude [°]".
    title : str, optional
        Plot title. Defaults to "Longitude Grid".
    marker : str, optional
        Marker style. Defaults to 'o'.
    **kwargs
        Additional keyword arguments passed to plot()

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : As input
        The matplotlib axes.
    """
    if fig is None:
        fig = plt.figure(figsize=(10, 8) if polar else (10, 6))
    
    if ax is None:
        ax = fig.add_subplot(1, 1, 1, polar=polar)
    
    if polar:
        # Convert longitude degrees to radians for polar plot
        # Longitude: -180 to 180 maps to full circle
        angles_rad = np.deg2rad(grid)
        radii = np.ones_like(grid)
        ax.plot(angles_rad, radii, marker=marker, **kwargs)
        ax.set_ylim(0, 1.2)
        ax.set_theta_zero_location("N")  # 0° at top (Prime Meridian)
        ax.set_theta_direction(-1)  # Clockwise (East positive)
        ax.set_title(title)
    else:
        indices = np.arange(len(grid))
        ax.plot(indices, grid, marker=marker, **kwargs)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_ylim(-180, 180)
        ax.grid(True, alpha=0.3)
        
        # Add horizontal lines at key longitudes
        ax.axhline(0, color='r', linestyle='--', alpha=0.3, label='Prime Meridian (0°)')
        ax.axhline(180, color='b', linestyle='--', alpha=0.3, label='Antimeridian (±180°)')
        ax.axhline(-180, color='b', linestyle='--', alpha=0.3)
        ax.legend()
    
    return fig, ax
