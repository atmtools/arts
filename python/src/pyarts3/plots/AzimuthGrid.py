""" Plotting routine for AzimuthGrid """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    grid: pyarts.arts.AzimuthGrid,
    *,
    fig=None,
    ax=None,
    polar: bool = False,
    xlabel: str = "Index",
    ylabel: str = "Azimuth Angle [°]",
    title: str = "Azimuth Grid",
    marker: str = 'o',
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

        pyarts.plots.AzimuthGrid.plot(azimuth, polar=True)

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
    xlabel : str, optional
        Label for x-axis (not used in polar). Defaults to "Index".
    ylabel : str, optional
        Label for y-axis (not used in polar). Defaults to "Azimuth Angle [°]".
    title : str, optional
        Plot title. Defaults to "Azimuth Grid".
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
        # Convert degrees to radians for polar plot
        # ARTS convention: 0° = North, 90° = East, 180° = South, 270° = West
        # Convert to math convention for polar plot: 0° = North at top
        angles_rad = np.deg2rad(grid)
        radii = np.ones_like(grid)
        ax.plot(angles_rad, radii, marker=marker, **kwargs)
        ax.set_ylim(0, 1.2)
        ax.set_theta_zero_location("N")  # 0° at North (top)
        ax.set_theta_direction(-1)  # Clockwise (East = 90° clockwise from North)
        ax.set_title(title)
    else:
        indices = np.arange(len(grid))
        ax.plot(indices, grid, marker=marker, **kwargs)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_ylim(0, 360)
        ax.grid(True, alpha=0.3)
        
        # Add horizontal lines at cardinal directions
        ax.axhline(0, color='r', linestyle='--', alpha=0.3, label='North (0°)')
        ax.axhline(90, color='g', linestyle='--', alpha=0.3, label='East (90°)')
        ax.axhline(180, color='b', linestyle='--', alpha=0.3, label='South (180°)')
        ax.axhline(270, color='m', linestyle='--', alpha=0.3, label='West (270°)')
        ax.legend()
    
    return fig, ax
