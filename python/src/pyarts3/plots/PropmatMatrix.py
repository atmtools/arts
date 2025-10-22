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
    component=None,
    element: tuple[int, int] | None = None,
    xlabel: str = "Column",
    ylabel: str = "Row",
    title: str = "Propagation Matrix",
    colorbar: bool = True,
    cmap: str = "viridis",
    **kwargs
):
    """
    Plot a propagation matrix as a 2D heatmap or grid of subplots.

    Modes:
    - Default (component=None, element=None): Grid of 16 subplots, one for each matrix element M[i,j].
    - Dot product (component=7-vector): Single heatmap of dot product with each cell.
    - Element mode (element=(i,j)): Single heatmap of M[i,j] over the domain.
    - Backward compatible: component=int plots compact component plane.

    Parameters
    ----------
    propmat_matrix : ~pyarts3.arts.PropmatMatrix
        A matrix of propagation matrices (7-vector or 4x4 each)
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes or list, optional
        Axes to plot on. If None, new axes are created.
    component : array-like or int or None, optional
        If None, show grid of 16 subplots (M[i,j]). If a 7-vector, plot dot product. If int, plot compact component plane.
    element : tuple[int, int] | None, optional
        If provided, plot heatmap of M[i,j] over the domain.
    xlabel, ylabel, title : str, optional
        Axis labels and plot title.
    colorbar : bool, optional
        Whether to show colorbar. Defaults to True.
    cmap : str, optional
        Colormap name. Defaults to "viridis".
    **kwargs
        Additional keyword arguments passed to imshow().

    Returns
    -------
    fig : Figure
        The matplotlib figure.
    ax : list or Axes
        List of axes (grid mode) or single axes (dot product/element/component mode).
    """
    if fig is None:
        fig = plt.figure(figsize=(10, 8))
    

    # Helper to expand a compact 7-vector [a,b,c,d,u,v,w] to 4x4 matrix
    def expand7(x):
        a, b, c, d, u, v, w = x
        return np.array([
            [a,  b,  c,  d],
            [b,  a,  u,  v],
            [c, -u,  a,  v],
            [d, -v, -w,  a],
        ])

    nx = len(propmat_matrix)
    ny = len(propmat_matrix[0]) if nx > 0 else 0

    if component is None and element is None:
        # Grid of plots: each subplot shows one matrix element (4x4)
        if ax is None:
            ax = []
            for i in range(4):
                row = []
                for j in range(4):
                    row.append(fig.add_subplot(4, 4, i*4 + j + 1))
                ax.append(row)
        for i in range(4):
            for j in range(4):
                data = np.empty((nx, ny))
                for x in range(nx):
                    for y in range(ny):
                        arr = np.asarray(propmat_matrix[x][y])
                        if arr.size == 7:
                            M = expand7(arr)
                            data[x, y] = M[i, j]
                        elif arr.shape == (4, 4):
                            data[x, y] = arr[i, j]
                        else:
                            data[x, y] = propmat_matrix[x][y][i, j]
                im = ax[i][j].imshow(data, aspect='auto', cmap=cmap, origin='lower', **kwargs)
                if colorbar:
                    plt.colorbar(im, ax=ax[i][j], label=f"M[{i},{j}]")
                ax[i][j].set_xlabel(xlabel)
                ax[i][j].set_ylabel(ylabel)
                ax[i][j].set_title(f"{title} M[{i},{j}]")
        return fig, ax
    elif component is not None and hasattr(component, '__len__') and len(component) == 7:
        # Dot product mode: component is a single 7-vector
        dot = np.empty((nx, ny))
        for x in range(nx):
            for y in range(ny):
                arr = np.asarray(propmat_matrix[x][y])
                if arr.size == 7:
                    dot[x, y] = np.dot(arr, component)
                elif arr.shape == (4, 4):
                    dot[x, y] = np.dot(arr.flatten(), expand7(component).flatten())
                else:
                    dot[x, y] = np.dot(np.asarray(propmat_matrix[x][y]).flatten(), np.asarray(component).flatten())
        if ax is None:
            ax = fig.add_subplot(1, 1, 1)
        im = ax.imshow(dot, aspect='auto', cmap=cmap, origin='lower', **kwargs)
        if colorbar:
            plt.colorbar(im, ax=ax, label="Dot product")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(f"{title} (dot product)")
        return fig, ax
    elif element is not None:
        i, j = element
        data = np.empty((nx, ny))
        for x in range(nx):
            for y in range(ny):
                arr = np.asarray(propmat_matrix[x][y])
                if arr.size == 7:
                    M = expand7(arr)
                    data[x, y] = M[i, j]
                elif arr.shape == (4, 4):
                    data[x, y] = arr[i, j]
                else:
                    data[x, y] = propmat_matrix[x][y][i, j]
        if ax is None:
            ax = fig.add_subplot(1, 1, 1)
        im = ax.imshow(data, aspect='auto', cmap=cmap, origin='lower', **kwargs)
        if colorbar:
            plt.colorbar(im, ax=ax, label=f"M[{i},{j}]")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(f"{title} M[{i},{j}]")
        return fig, ax
    else:
        # Back-compat: plot compact component plane
        if component is None:
            component = 0
        data = np.asarray(propmat_matrix[:, :, component])
        if ax is None:
            ax = fig.add_subplot(1, 1, 1)
        im = ax.imshow(data, aspect='auto', cmap=cmap, origin='lower', **kwargs)
        if colorbar:
            plt.colorbar(im, ax=ax, label=f"Component {component}")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(f"{title} (Component {component})")
        return fig, ax
