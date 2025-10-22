""" Plotting routine for AscendingGrid """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    grid: pyarts.arts.AscendingGrid,
    *,
    fig=None,
    ax=None,
    xlabel: str = "Index",
    ylabel: str = "Grid Value",
    title: str = "Ascending Grid",
    marker: str = 'o',
    show_info: bool = True,
    **kwargs
):
    """Plot an AscendingGrid as a scatter/line plot.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create a frequency grid
        freqs = pyarts.arts.AscendingGrid(np.logspace(9, 12, 20))

        pyarts.plots.AscendingGrid.plot(freqs, ylabel="Frequency [Hz]", 
                                        title="Frequency Grid")

    Parameters
    ----------
    grid : ~pyarts3.arts.AscendingGrid
        A sorted ascending grid of values
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    xlabel : str, optional
        Label for x-axis. Defaults to "Index".
    ylabel : str, optional
        Label for y-axis. Defaults to "Grid Value".
    title : str, optional
        Plot title. Defaults to "Ascending Grid".
    marker : str, optional
        Marker style. Defaults to 'o'.
    show_info : bool, optional
        Whether to display info text box with statistics. Defaults to True.
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
        fig = plt.figure(figsize=(10, 6))
    
    if ax is None:
        ax = fig.add_subplot(1, 1, 1)
    
    indices = np.arange(len(grid))
    ax.plot(indices, grid, marker=marker, **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    
    # Add statistics as text
    if show_info:
        info_text = f"Points: {len(grid)}\n"
        info_text += f"Min: {np.min(grid):.2e}\n"
        info_text += f"Max: {np.max(grid):.2e}\n"
        info_text += f"Range: {np.max(grid) - np.min(grid):.2e}"
        
        ax.text(0.02, 0.98, info_text, transform=ax.transAxes,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    return fig, ax
