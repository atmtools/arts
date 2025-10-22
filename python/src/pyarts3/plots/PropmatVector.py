""" Plotting routine for PropmatVector """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

__all__ = [
    'plot',
]


def plot(
    propmat_vector: pyarts.arts.PropmatVector,
    *,
    fig=None,
    ax=None,
    freqs: np.ndarray | None = None,
    component=None,
    element: tuple[int, int] | None = None,
    xlabel: str = "Index",
    ylabel: str = "Propagation Matrix Component",
    title: str = "Propagation Matrix Vector",
    **kwargs
):
    """
    Plot a propagation matrix vector.

    Modes:
    - Default (component=None, element=None): Grid of 16 subplots, one for each matrix element M[i,j].
    - Dot product (component=7-vector): Single plot of dot product with each sample.
    - Element mode (element=(i,j)): Single plot of M[i,j] across the vector.
    - Backward compatible: component=int plots compact component.

    Parameters
    ----------
    propmat_vector : ~pyarts3.arts.PropmatVector
        A vector of propagation matrices (7-vector or 4x4 each)
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes or list, optional
        Axes to plot on. If None, new axes are created.
    freqs : :class:`~numpy.ndarray` | None, optional
        Frequency or position grid for x-axis. If None, uses indices.
    component : array-like or int or None, optional
        If None, show grid of 16 subplots (M[i,j]). If a 7-vector, plot dot product. If int, plot compact component.
    element : tuple[int, int] | None, optional
        If provided, plot M[i,j] across the vector.
    xlabel, ylabel, title : str, optional
        Axis labels and plot title.
    **kwargs
        Additional keyword arguments passed to matplotlib plot().

    Returns
    -------
    fig : Figure
        The matplotlib figure.
    ax : list or Axes
        List of axes (grid mode) or single axes (dot product/element/component mode).
    """
    if fig is None:
        fig = plt.figure(figsize=(10, 6))
    
    if ax is None:
        ax = fig.add_subplot(1, 1, 1)
    
    # Helper to expand a compact 7-vector [a,b,c,d,u,v,w] to 4x4 matrix
    def expand7(x):
        a, b, c, d, u, v, w = x
        return np.array([
            [a,  b,  c,  d],
            [b,  a,  u,  v],
            [c, -u,  a,  v],
            [d, -v, -w,  a],
        ])


    n = len(propmat_vector)

    if component is None and element is None:
        # Grid of plots: each subplot shows one matrix element (4x4)
        # Normalize axes to a 4x4 grid (create if needed)
        def to_grid_axes(ax_in):
            if ax_in is None:
                grid = []
                for i in range(4):
                    row = []
                    for j in range(4):
                        row.append(fig.add_subplot(4, 4, i*4 + j + 1))
                    grid.append(row)
                return grid
            if isinstance(ax_in, Axes):
                # Provided a single axes; create a full grid for grid mode
                return to_grid_axes(None)
            # List or array-like
            try:
                # nested list-of-lists
                if len(ax_in) == 4 and all(hasattr(ax_in[i], '__len__') and len(ax_in[i]) == 4 for i in range(4)):
                    return ax_in
            except Exception:
                pass
            try:
                # flat list of 16 axes
                if hasattr(ax_in, '__len__') and len(ax_in) == 16:
                    return [list(ax_in[i*4:(i+1)*4]) for i in range(4)]
            except Exception:
                pass
            # Fallback: create new grid
            return to_grid_axes(None)

        ax = to_grid_axes(ax)
        for i in range(4):
            for j in range(4):
                series = []
                for m in propmat_vector:
                    arr = np.asarray(m)
                    if arr.size == 7:
                        M = expand7(arr)
                        series.append(M[i, j])
                    elif arr.shape == (4, 4):
                        series.append(arr[i, j])
                    else:
                        series.append(m[i, j])
                data = np.asarray(series)
                f = freqs if freqs is not None else np.arange(n)
                ax[i][j].plot(f, data, **kwargs)
                ax[i][j].set_xlabel(xlabel)
                ax[i][j].set_ylabel(ylabel)
                ax[i][j].set_title(f"{title} M[{i},{j}]")
                ax[i][j].grid(True, alpha=0.3)
        return fig, ax
    elif component is not None and hasattr(component, '__len__') and len(component) == 7:
        # Dot product mode: component is a single 7-vector
        dot = []
        for m in propmat_vector:
            arr = np.asarray(m)
            if arr.size == 7:
                dot.append(np.dot(arr, component))
            elif arr.shape == (4, 4):
                dot.append(np.dot(arr.flatten(), expand7(component).flatten()))
            else:
                dot.append(np.dot(np.asarray(m).flatten(), np.asarray(component).flatten()))
        data = np.asarray(dot)
        f = freqs if freqs is not None else np.arange(n)
        # Normalize to a single axes
        if ax is None or (hasattr(ax, '__len__') and len(ax) != 0):
            # If ax is a list (grid from previous call), pick the first
            try:
                ax = ax[0][0] if hasattr(ax[0], '__len__') else ax[0]
            except Exception:
                ax = fig.add_subplot(1, 1, 1)
        elif not isinstance(ax, Axes):
            ax = fig.add_subplot(1, 1, 1)
        ax.plot(f, data, **kwargs)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(f"{title} (dot product)")
        ax.grid(True, alpha=0.3)
        return fig, ax
    elif element is not None:
        # Plot a particular matrix element (i,j) across vector
        i, j = element
        series = []
        for m in propmat_vector:
            arr = np.asarray(m)
            if arr.size == 7:
                M = expand7(arr)
                series.append(M[i, j])
            elif arr.shape == (4, 4):
                series.append(arr[i, j])
            else:
                series.append(m[i, j])
        data = np.asarray(series)
        f = freqs if freqs is not None else np.arange(n)
        # Normalize to a single axes
        if ax is None or (hasattr(ax, '__len__') and len(ax) != 0):
            try:
                ax = ax[0][0] if hasattr(ax[0], '__len__') else ax[0]
            except Exception:
                ax = fig.add_subplot(1, 1, 1)
        elif not isinstance(ax, Axes):
            ax = fig.add_subplot(1, 1, 1)
        ax.plot(f, data, **kwargs)
        ax.set_title(f"{title} (M[{i},{j}])")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)
        return fig, ax
    else:
        # Back-compat: plot a single compact component (0-6)
        if component is None:
            component = 0
        data = np.asarray([np.asarray(m)[component] for m in propmat_vector])
        f = freqs if freqs is not None else np.arange(n)
        if ax is None:
            ax = fig.add_subplot(1, 1, 1)
        ax.plot(f, data, **kwargs)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(f"{title} (Component {component})")
        ax.grid(True, alpha=0.3)
        return fig, ax
