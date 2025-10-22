""" Plotting routine for LatGrid """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    grid: pyarts.arts.LatGrid,
    *,
    fig=None,
    ax=None,
    polar: bool = False,
    xlabel: str = "Index",
    ylabel: str = "Latitude [°]",
    title: str = "Latitude Grid",
    marker: str = 'o',
    **kwargs
):
    """Plot a LatGrid showing latitude values.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create a latitude grid
        lats = pyarts.arts.LatGrid(np.linspace(-90, 90, 19))

        pyarts.plots.LatGrid.plot(lats, polar=True)

    Parameters
    ----------
    grid : ~pyarts3.arts.LatGrid
        A sorted grid of latitude values [-90, 90]
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    polar : bool, optional
        If True, use polar plot. Defaults to False.
    xlabel : str, optional
        Label for x-axis (not used in polar). Defaults to "Index".
    ylabel : str, optional
        Label for y-axis (not used in polar). Defaults to "Latitude [°]".
    title : str, optional
        Plot title. Defaults to "Latitude Grid".
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
        # Convert latitudes to polar coordinates
        # Display: 90° (North Pole) at top, -90° (South Pole) at bottom
        # Need to flip the angle mapping: latitude → (90 - latitude)
        # So: lat=90° → angle=0° (top), lat=0° → angle=90° (right), lat=-90° → angle=180° (bottom)
        flipped_angles = 90 - grid
        angles_rad = np.deg2rad(flipped_angles)
        radii = np.ones_like(grid)
        ax.plot(angles_rad, radii, marker=marker, **kwargs)
        ax.set_ylim(0, 1.2)
        ax.set_theta_zero_location("N")  # 0° at top (North Pole, lat=90°)
        ax.set_theta_direction(1)  # Counter-clockwise
        ax.set_thetamin(0)  # Limit to half circle
        ax.set_thetamax(180)
        
        # Relabel the angles to show latitude values
        theta_labels = ['90°N', '45°N', '0°', '45°S', '90°S']
        ax.set_thetagrids([0, 45, 90, 135, 180], labels=theta_labels)
        ax.set_title(title)
    else:
        indices = np.arange(len(grid))
        ax.plot(indices, grid, marker=marker, **kwargs)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_ylim(-90, 90)
        ax.grid(True, alpha=0.3)
        
        # Add horizontal lines at key latitudes
        ax.axhline(0, color='r', linestyle='--', alpha=0.3, label='Equator (0°)')
        ax.axhline(90, color='b', linestyle='--', alpha=0.3, label='North Pole (90°)')
        ax.axhline(-90, color='b', linestyle='--', alpha=0.3, label='South Pole (-90°)')
        ax.legend()
    
    return fig, ax
