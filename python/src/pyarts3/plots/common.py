import matplotlib.pyplot as plt
import numpy as np


def default_fig_ax(fig=None, ax=None, nrows=1, ncols=1, N=-1, fig_kwargs={}, ax_kwargs={}):
    """Utility to create default matplotlib figure and axes if not provided.

    Parameters
    ----------
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    fig_kwargs : dict, optional
        Keyword arguments for creating new figure if fig is None. Defaults to {}.
    ax_kwargs : dict, optional
        Keyword arguments for creating new axes if ax is None. Defaults to {}.

    Returns
    -------
    fig : Figure
        The matplotlib figure.
    ax : Axes
        The matplotlib axes.
    """
    fig = plt.figure(**fig_kwargs) if fig is None else fig
    if ax is None:
        ax = fig.subplots(nrows, ncols, **ax_kwargs)
        if N >= 0:
            for a in ax.flatten()[N:]:
                fig.delaxes(a)

    return fig, ax


def select_flat_ax(ax, index):
    """Utility to select a single Axes from possibly multi-dimensional axes array.

    Parameters
    ----------
    ax : Axes or array of Axes
        The matplotlib axes or array of axes.
    index : int
        The index of the desired Axes in flattened order.

    Returns
    -------
    selected_ax : Axes
        The selected matplotlib axes.
    """
    if isinstance(ax, plt.Axes):
        return ax
    elif isinstance(ax, np.ndarray):
        return ax.flat[index]
    elif isinstance(ax, (list, tuple)):
        return ax[index]
    else:
        return ax.flatten()[index]
