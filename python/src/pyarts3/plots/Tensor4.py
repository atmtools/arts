""" Plotting routine for Tensor4 """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    tensor: pyarts.arts.Tensor4,
    *,
    fig=None,
    slice_dims: tuple = (0, 1),
    slice_indices: tuple = (0, 0),
    xlabel: str = "Dimension 3",
    ylabel: str = "Dimension 2",
    title: str | None = None,
    colorbar: bool = True,
    cmap: str = "viridis",
    **kwargs
):
    """Plot a 2D slice of a Tensor4 as a heatmap.

    Parameters
    ----------
    tensor : ~pyarts3.arts.Tensor4
        A 4D array of numeric values
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    slice_dims : tuple, optional
        Which two dimensions to slice (e.g., (0, 1)). Defaults to (0, 1).
    slice_indices : tuple, optional
        Indices for the sliced dimensions. Defaults to (0, 0).
    xlabel : str, optional
        Label for x-axis. Defaults to "Dimension 3".
    ylabel : str, optional
        Label for y-axis. Defaults to "Dimension 2".
    title : str | None, optional
        Plot title. If None, auto-generates. Defaults to None.
    colorbar : bool, optional
        Whether to show colorbar. Defaults to True.
    cmap : str, optional
        Colormap name. Defaults to "viridis".
    **kwargs
        Additional keyword arguments passed to imshow()

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : As input
        The matplotlib axes.
    """
    if fig is None:
        fig = plt.figure(figsize=(10, 8))
    
    ax = fig.add_subplot(1, 1, 1)
    
    # Extract 2D slice by fixing two dimensions
    shape = tensor.shape
    idx = [slice(None)] * 4
    for i, (dim, val) in enumerate(zip(slice_dims, slice_indices)):
        idx[dim] = val
    
    data = tensor[tuple(idx)]
    
    # Squeeze to 2D if needed
    while data.ndim > 2:
        data = data[0]
    
    im = ax.imshow(data, aspect='auto', cmap=cmap, origin='lower', **kwargs)
    
    if colorbar:
        plt.colorbar(im, ax=ax)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    if title is None:
        title = f"Tensor4 Slice (dims={slice_dims}, indices={slice_indices})"
    ax.set_title(title)
    
    return fig, ax
