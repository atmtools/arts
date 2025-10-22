""" Plotting routine for SortedGriddedField3 """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    gridded_field: pyarts.arts.SortedGriddedField3,
    *,
    fig=None,
    ax=None,
    slice_dim: int = 0,
    slice_idx: int = 0,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
    colorbar: bool = True,
    cmap: str = "viridis",
    **kwargs
):
    """Plot a 2D slice of SortedGriddedField3 as a heatmap.

    Parameters
    ----------
    gridded_field : ~pyarts3.arts.SortedGriddedField3
        A 3D sorted gridded field with ascending grids
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    slice_dim : int, optional
        Which dimension to slice (0, 1, or 2). Defaults to 0.
    slice_idx : int, optional
        Index along the slice dimension. Defaults to 0.
    xlabel : str | None, optional
        Label for x-axis. If None, uses grid name. Defaults to None.
    ylabel : str | None, optional
        Label for y-axis. If None, uses grid name. Defaults to None.
    title : str | None, optional
        Plot title. If None, auto-generates. Defaults to None.
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
    
    # Extract 2D slice
    if slice_dim == 0:
        data = gridded_field.data[slice_idx, :, :]
        xgrid = gridded_field.grids[1]
        ygrid = gridded_field.grids[2]
        xname = gridded_field.gridnames[1]
        yname = gridded_field.gridnames[2]
        slice_value = gridded_field.grids[0][slice_idx]
        slice_name = gridded_field.gridnames[0]
    elif slice_dim == 1:
        data = gridded_field.data[:, slice_idx, :]
        xgrid = gridded_field.grids[0]
        ygrid = gridded_field.grids[2]
        xname = gridded_field.gridnames[0]
        yname = gridded_field.gridnames[2]
        slice_value = gridded_field.grids[1][slice_idx]
        slice_name = gridded_field.gridnames[1]
    else:  # slice_dim == 2
        data = gridded_field.data[:, :, slice_idx]
        xgrid = gridded_field.grids[0]
        ygrid = gridded_field.grids[1]
        xname = gridded_field.gridnames[0]
        yname = gridded_field.gridnames[1]
        slice_value = gridded_field.grids[2][slice_idx]
        slice_name = gridded_field.gridnames[2]
    
    im = ax.pcolormesh(xgrid, ygrid, data, cmap=cmap, shading='auto', **kwargs)
    
    if colorbar:
        plt.colorbar(im, ax=ax, label=gridded_field.dataname)
    
    ax.set_xlabel(xlabel if xlabel is not None else xname)
    ax.set_ylabel(ylabel if ylabel is not None else yname)
    
    if title is None:
        title = f"{gridded_field.dataname} ({slice_name}={slice_value:.2e})"
    ax.set_title(title)
    
    return fig, ax
