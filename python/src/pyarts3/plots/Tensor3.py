""" Plotting routine for Tensor3 """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    tensor: pyarts.arts.Tensor3,
    *,
    fig=None,
    slice_dim: int = 0,
    slice_idx: int = 0,
    xlabel: str = "Dimension 2",
    ylabel: str = "Dimension 1",
    title: str | None = None,
    colorbar: bool = True,
    cmap: str = "viridis",
    **kwargs
):
    """Plot a slice of a Tensor3 as a 2D heatmap.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create a 3D tensor
        x = np.linspace(-2, 2, 20)
        y = np.linspace(-2, 2, 25)
        z = np.linspace(-2, 2, 15)
        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
        data = np.exp(-(X**2 + Y**2 + Z**2))
        
        tensor = pyarts.arts.Tensor3(data)

        pyarts.plots.Tensor3.plot(tensor, slice_dim=2, slice_idx=7)

    Parameters
    ----------
    tensor : ~pyarts3.arts.Tensor3
        A 3D array of numeric values
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    slice_dim : int, optional
        Which dimension to slice (0, 1, or 2). Defaults to 0.
    slice_idx : int, optional
        Index along the slice dimension. Defaults to 0.
    xlabel : str, optional
        Label for x-axis. Defaults to "Dimension 2".
    ylabel : str, optional
        Label for y-axis. Defaults to "Dimension 1".
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
    
    # Extract 2D slice
    if slice_dim == 0:
        data = tensor[slice_idx, :, :]
    elif slice_dim == 1:
        data = tensor[:, slice_idx, :]
    else:  # slice_dim == 2
        data = tensor[:, :, slice_idx]
    
    im = ax.imshow(data, aspect='auto', cmap=cmap, origin='lower', **kwargs)
    
    if colorbar:
        plt.colorbar(im, ax=ax)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    if title is None:
        title = f"Tensor3 Slice (dim={slice_dim}, idx={slice_idx})"
    ax.set_title(title)
    
    return fig, ax
