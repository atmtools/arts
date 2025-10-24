""" Plotting routine for ZenithGrid """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    grid: pyarts.arts.ZenithGrid,
    *,
    fig=None,
    ax=None,
    polar: bool = False,
    xlabel: str = "Index",
    ylabel: str = "Zenith Angle [°]",
    title: str = "Zenith Grid",
    marker: str = 'o',
    **kwargs
):
    """Plot a ZenithGrid showing zenith angles.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create zenith angles from 0° (up) to 180° (down)
        zenith = pyarts.arts.ZenithGrid(np.linspace(0, 180, 19))

        pyarts.plots.ZenithGrid.plot(zenith, polar=True)

    Parameters
    ----------
    grid : ~pyarts3.arts.ZenithGrid
        A sorted grid of zenith angles [0, 180]
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    polar : bool, optional
        If True, use polar plot. Defaults to False.
    xlabel : str, optional
        Label for x-axis (not used in polar). Defaults to "Index".
    ylabel : str, optional
        Label for y-axis (not used in polar). Defaults to "Zenith Angle [°]".
    title : str, optional
        Plot title. Defaults to "Zenith Grid".
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
        # Zenith 0° = up, 180° = down (half circle)
        angles_rad = np.deg2rad(grid)
        radii = np.ones_like(grid)
        ax.plot(angles_rad, radii, marker=marker, **kwargs)
        ax.set_ylim(0, 1.2)
        ax.set_theta_zero_location("N")  # 0° at top
        ax.set_theta_direction(-1)  # Clockwise
        ax.set_thetamin(0)  # Limit to half circle
        ax.set_thetamax(180)
        ax.set_title(title)
    else:
        indices = np.arange(len(grid))
        ax.plot(indices, grid, marker=marker, **kwargs)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_ylim(0, 180)
        ax.grid(True, alpha=0.3)
        
        # Add horizontal lines at key angles
        ax.axhline(0, color='r', linestyle='--', alpha=0.3, label='Zenith (0°)')
        ax.axhline(90, color='g', linestyle='--', alpha=0.3, label='Horizon (90°)')
        ax.axhline(180, color='b', linestyle='--', alpha=0.3, label='Nadir (180°)')
        ax.legend()
    
    return fig, ax
