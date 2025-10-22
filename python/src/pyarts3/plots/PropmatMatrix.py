""" Plotting routine for PropmatMatrix """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    propmat_matrix: pyarts.arts.PropmatMatrix,
    *,
    fig=None,
    ax=None,
    component: int = 0,
    xlabel: str = "Column",
    ylabel: str = "Row",
    title: str = "Propagation Matrix",
    colorbar: bool = True,
    cmap: str = "viridis",
    **kwargs
):
    """Plot a propagation matrix as a 2D heatmap for a specific component.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create a 2D absorption field
        nx, ny = 30, 25
        x = np.linspace(-2, 2, nx)
        y = np.linspace(-2, 2, ny)
        X, Y = np.meshgrid(x, y, indexing='ij')
        
        propmat_mat = pyarts.arts.PropmatMatrix(np.zeros((nx, ny, 7)))
        propmat_mat[:, :, 0] = np.exp(-(X**2 + Y**2))  # Absorption component

        pyarts.plots.PropmatMatrix.plot(propmat_mat, component=0)

    Parameters
    ----------
    propmat_matrix : ~pyarts3.arts.PropmatMatrix
        A matrix of propagation matrices (7 components each)
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    component : int, optional
        Which component to plot (0-6). Defaults to 0.
    xlabel : str, optional
        Label for x-axis. Defaults to "Column".
    ylabel : str, optional
        Label for y-axis. Defaults to "Row".
    title : str, optional
        Plot title. Defaults to "Propagation Matrix".
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
    
    if ax is None:
        ax = fig.add_subplot(1, 1, 1)
    
    # Extract component from propmat matrix
    data = propmat_matrix[:, :, component]
    
    im = ax.imshow(data, aspect='auto', cmap=cmap, origin='lower', **kwargs)
    
    if colorbar:
        plt.colorbar(im, ax=ax, label=f"Component {component}")
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"{title} (Component {component})")
    
    return fig, ax
