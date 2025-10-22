""" Plotting routine for StokvecVector """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes

__all__ = [
    'plot',
]


def plot(
    stokvec_vector: pyarts.arts.StokvecVector,
    *,
    fig=None,
    ax=None,
    freqs: np.ndarray | None = None,
    component=None,
    xlabel: str = "Index",
    ylabel: str = "Value",
    title: str = "Stokes Vector Components",
    labels: list[str] = ['I', 'Q', 'U', 'V'],
    **kwargs
):
    """
    Plot a vector of Stokes vectors.

    Modes:
    - Default (component=None): Grid of 4 subplots, one for each Stokes component (I, Q, U, V).
    - Dot product (component=4-vector): Single plot of dot product with each sample.

    Parameters
    ----------
    stokvec_vector : ~pyarts3.arts.StokvecVector
        A vector of Stokes vectors (each with 4 components: I, Q, U, V)
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes or list, optional
        Axes to plot on. If None, new axes are created.
    freqs : :class:`~numpy.ndarray` | None, optional
        Frequency or position grid for x-axis. If None, uses indices.
    component : array-like or None, optional
        If None, show grid of 4 subplots (I,Q,U,V). If a 4-vector, plot dot product with each sample.
    xlabel, ylabel, title, labels : str or list, optional
        Axis labels and plot title.
    **kwargs
        Additional keyword arguments passed to matplotlib plot().

    Returns
    -------
    fig : Figure
        The matplotlib figure.
    ax : list or Axes
        List of axes (grid mode) or single axes (dot product mode).
    """
    if fig is None:
        fig = plt.figure(figsize=(12, 10), constrained_layout=True)
    
    if ax is None:
        ax = []
        for i in range(4):
            ax.append(fig.add_subplot(2, 2, i + 1))
    
    if freqs is None:
        freqs = np.arange(len(stokvec_vector))

    n = len(stokvec_vector)

    if component is None:
        # Grid of plots: each subplot shows one Stokes component
        def to_grid_axes(ax_in):
            if ax_in is None:
                return [fig.add_subplot(2, 2, i + 1) for i in range(4)]
            if isinstance(ax_in, Axes):
                return to_grid_axes(None)
            try:
                if hasattr(ax_in, '__len__') and len(ax_in) == 4 and not hasattr(ax_in[0], '__len__'):
                    return ax_in
            except Exception:
                pass
            return to_grid_axes(None)
        ax = to_grid_axes(ax)
        # Robust extraction: support wrapped 4x1 objects
        data_ok = False
        try:
            arr = np.asarray(stokvec_vector)
            if arr.ndim == 2 and arr.shape[1] == 4 and len(arr) == n:
                comps = [arr[:, i] for i in range(4)]
                data_ok = True
            else:
                data_ok = False
        except Exception:
            data_ok = False
        if not data_ok:
            comps = [np.array([sv[i] for sv in stokvec_vector]) for i in range(4)]
        f = freqs if freqs is not None else np.arange(n)
        for i, label in enumerate(labels):
            subplot = ax[i] if isinstance(ax, list) else ax
            data = comps[i]
            subplot.plot(f, data, label=label, **kwargs)
            subplot.set_xlabel(xlabel)
            subplot.set_ylabel(ylabel)
            subplot.set_title(f"{title} - {label}")
            subplot.legend()
            subplot.grid(True, alpha=0.3)
        return fig, ax
    elif hasattr(component, '__len__') and len(component) == 4:
        # Dot product mode: component is a single 4-vector
        dot = []
        for sv in stokvec_vector:
            arr = np.asarray(sv)
            if arr.size == 4:
                dot.append(np.dot(arr, component))
            else:
                dot.append(np.dot(np.array([sv[i] for i in range(4)]), component))
        data = np.asarray(dot)
        f = freqs if freqs is not None else np.arange(n)
        # Normalize to a single axes (pick first if a list of axes was provided)
        if ax is None or (hasattr(ax, '__len__') and len(ax) != 0):
            try:
                ax = ax[0]
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
