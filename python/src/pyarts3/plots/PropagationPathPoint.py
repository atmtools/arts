""" Plotting routine for PropagationPathPoint """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    path_point: pyarts.arts.PropagationPathPoint,
    *,
    fig=None,
    ax=None,
    show_info: bool = True,
    **kwargs,
):
    """Plot a propagation path point showing position on a map.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts

        # Create a path point
        path = pyarts.arts.PropagationPathPoint()
        path.pos = [100e3, 45.0, 10.0]  # altitude, lat, lon
        path.los = [120.0, 45.0]  # zenith, azimuth

        pyarts.plots.PropagationPathPoint.plot(path)

    Parameters
    ----------
    path_point : ~pyarts3.arts.PropagationPathPoint
        A propagation path point with position and line-of-sight
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    show_info : bool, optional
        Whether to display info text box with position and LOS details. Defaults to True.
    **kwargs
        Additional keyword arguments passed to matplotlib functions (ax.plot, ax.arrow, etc.).

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : As input
        The matplotlib axes.
    """
    if fig is None:
        fig = plt.figure(figsize=(10, 8))
    
    if ax is None:
        ax = fig.add_subplot(1, 1, 1)
    
    # Extract position: [altitude, latitude, longitude]
    pos = path_point.pos
    altitude = pos[0]  # meters
    latitude = pos[1]   # degrees
    longitude = pos[2]  # degrees
    
    # Extract line-of-sight: [zenith, azimuth]
    # Azimuth convention: 0° = North, 90° = East, 180° = South, 270° = West
    los = path_point.los
    zenith = los[0]    # degrees
    azimuth = los[1]   # degrees (0° = North, 90° = East)
    
    # Plot position on map with a big cross
    ax.plot(longitude, latitude, 'rx', markersize=20, markeredgewidth=3, label='Position', **kwargs)
    
    # Add arrow showing line-of-sight direction (horizontal projection)
    # Convert azimuth to standard mathematical angle: 0° = North → 90° math angle
    # ARTS: 0°=N, 90°=E, 180°=S, 270°=W
    # Math: 0°=E, 90°=N, 180°=W, 270°=S
    # Conversion: math_angle = 90 - azimuth
    arrow_length = 2.0  # degrees
    math_angle_rad = np.deg2rad(90 - azimuth)
    dx = arrow_length * np.cos(math_angle_rad)
    dy = arrow_length * np.sin(math_angle_rad)
    
    ax.arrow(longitude, latitude, dx, dy, 
             head_width=0.5, head_length=0.3, fc='blue', ec='blue', 
             linewidth=2, alpha=0.7, label='LOS direction', **kwargs)
    
    # Set reasonable axis limits around the point
    lat_range = 5.0  # degrees
    lon_range = 5.0  # degrees
    ax.set_xlim(longitude - lon_range, longitude + lon_range)
    ax.set_ylim(latitude - lat_range, latitude + lat_range)
    
    ax.set_xlabel('Longitude [°]')
    ax.set_ylabel('Latitude [°]')
    ax.set_title('Propagation Path Point Position')
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal', adjustable='box')
    ax.legend()
    
    # Add info text with position and line-of-sight details
    if show_info:
        # Add cardinal direction for azimuth
        if azimuth < 22.5 or azimuth >= 337.5:
            direction = "N"
        elif azimuth < 67.5:
            direction = "NE"
        elif azimuth < 112.5:
            direction = "E"
        elif azimuth < 157.5:
            direction = "SE"
        elif azimuth < 202.5:
            direction = "S"
        elif azimuth < 247.5:
            direction = "SW"
        elif azimuth < 292.5:
            direction = "W"
        else:
            direction = "NW"
        
        info_text = f"Altitude: {altitude:.2e} m\n"
        info_text += f"Latitude: {latitude:.4f}°\n"
        info_text += f"Longitude: {longitude:.4f}°\n"
        info_text += f"Zenith: {zenith:.2f}°\n"
        info_text += f"Azimuth: {azimuth:.2f}° ({direction})"
        
        ax.text(0.02, 0.98, info_text, transform=ax.transAxes,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    return fig, ax
