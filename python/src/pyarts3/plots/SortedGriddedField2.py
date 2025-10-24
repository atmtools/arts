""" Plotting routine for SortedGriddedField2 """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    gridded_field: pyarts.arts.SortedGriddedField2,
    *,
    fig=None,
    ax=None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
    colorbar: bool = True,
    cmap: str = "viridis",
    **kwargs
):
    """Plot a SortedGriddedField2 as a 2D heatmap using its sorted grids.

    Parameters
    ----------
    gridded_field : ~pyarts3.arts.SortedGriddedField2
        A 2D sorted gridded field with ascending grids
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    xlabel : str | None, optional
        Label for x-axis. If None, uses grid name. Defaults to None.
    ylabel : str | None, optional
        Label for y-axis. If None, uses grid name. Defaults to None.
    title : str | None, optional
        Plot title. If None, uses data name. Defaults to None.
    colorbar : bool, optional
        Whether to show colorbar. Defaults to True.
    cmap : str, optional
        Colormap name. Defaults to "viridis".
    **kwargs
        Additional keyword arguments passed to pcolormesh()

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
    
    # Get grids and data
    xgrid = gridded_field.grids[0]
    ygrid = gridded_field.grids[1]
    data = gridded_field.data
    
    # Plot using pcolormesh
    im = ax.pcolormesh(xgrid, ygrid, data, cmap=cmap, **kwargs)
    
    if colorbar:
        plt.colorbar(im, ax=ax, label=gridded_field.dataname if title is None else title)
    
    ax.set_xlabel(xlabel if xlabel is not None else gridded_field.gridnames[0])
    ax.set_ylabel(ylabel if ylabel is not None else gridded_field.gridnames[1])
    ax.set_title(title if title is not None else gridded_field.dataname)
    
    return fig, ax
