""" Plotting routine for StokvecMatrix """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    stokvec_matrix: pyarts.arts.StokvecMatrix,
    *,
    fig=None,
    ax=None,
    component: int = 0,
    xlabel: str = "Column",
    ylabel: str = "Row",
    title: str = "Stokes Vector Matrix",
    colorbar: bool = True,
    cmap: str = "viridis",
    labels: list[str] = ['I', 'Q', 'U', 'V'],
    **kwargs
):
    """Plot a Stokes vector matrix as a 2D heatmap for a specific component.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create a 2D radiance field
        nx, ny = 40, 30
        x = np.linspace(0, 2*np.pi, nx)
        y = np.linspace(0, 2*np.pi, ny)
        X, Y = np.meshgrid(x, y, indexing='ij')
        
        stokes_mat = pyarts.arts.StokvecMatrix(np.zeros((nx, ny, 4)))
        stokes_mat[:, :, 0] = np.sin(X) * np.cos(Y)  # I component

        pyarts.plots.StokvecMatrix.plot(stokes_mat, component=0)

    Parameters
    ----------
    stokvec_matrix : ~pyarts3.arts.StokvecMatrix
        A matrix of Stokes vectors (4 components each)
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    component : int, optional
        Which Stokes component to plot (0=I, 1=Q, 2=U, 3=V). Defaults to 0.
    xlabel : str, optional
        Label for x-axis. Defaults to "Column".
    ylabel : str, optional
        Label for y-axis. Defaults to "Row".
    title : str, optional
        Plot title. Defaults to "Stokes Vector Matrix".
    colorbar : bool, optional
        Whether to show colorbar. Defaults to True.
    cmap : str, optional
        Colormap name. Defaults to "viridis".
    labels : list[str], optional
        Labels for the 4 Stokes parameters. Defaults to ['I', 'Q', 'U', 'V'].
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
    
    # Extract component from stokvec matrix
    data = stokvec_matrix[:, :, component]
    
    im = ax.imshow(data, aspect='auto', cmap=cmap, origin='lower', **kwargs)
    
    if colorbar:
        plt.colorbar(im, ax=ax, label=f"{labels[component]} component")
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"{title} - {labels[component]}")
    
    return fig, ax
