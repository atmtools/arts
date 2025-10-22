""" Plotting routine for Time """

import pyarts3 as pyarts
import matplotlib.pyplot as plt
from datetime import datetime

__all__ = [
    'plot',
]


def plot(
    time: pyarts.arts.Time,
    *,
    fig=None,
    ax=None,
    title: str = "Time",
    format: str = "%Y-%m-%d %H:%M:%S",
    show_info: bool = True,
    **kwargs,
):
    """Display a Time object as a formatted text plot.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts

        # Create a time object (current date for example)
        time = pyarts.arts.Time("2025-10-22 12:30:45")

        pyarts.plots.Time.plot(time)

    Parameters
    ----------
    time : ~pyarts3.arts.Time
        A time stamp
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    title : str, optional
        Plot title. Defaults to "Time".
    format : str, optional
        Time format string. Defaults to "%Y-%m-%d %H:%M:%S".
    show_info : bool, optional
        Whether to display the time info text box. Defaults to True.
    **kwargs
        Additional keyword arguments passed to matplotlib text function.

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : As input
        The matplotlib axes.
    """
    if fig is None:
        fig = plt.figure(figsize=(8, 4))
    
    if ax is None:
        ax = fig.add_subplot(1, 1, 1)
    
    # Hide axes
    ax.axis('off')
    
    if show_info:
        # Convert time to string representation
        try:
            # Try to convert to datetime for formatting
            time_str = str(time)
        except:
            time_str = repr(time)
        
        # Display as large centered text
        ax.text(0.5, 0.5, time_str, 
                transform=ax.transAxes,
                fontsize=20,
                verticalalignment='center',
                horizontalalignment='center',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8),
                **kwargs)
    
    ax.set_title(title)
    
    return fig, ax
