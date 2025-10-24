""" Plotting routine for AtmPoint """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    atm_point: pyarts.arts.AtmPoint,
    *,
    fig=None,
    ax=None,
    keys: list[str] | None = None,
    **kwargs
):
    """Plot atmospheric point parameters as a bar chart.
    
    Note: This function creates multiple subplots, so the ax parameter is accepted
    for API consistency but not used.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts

        ws = pyarts.Workspace()
        ws.atmospheric_fieldRead(toa=100e3, basename="planets/Earth/afgl/tropical/")
        
        # Extract atmospheric point at 50 km altitude
        atm_point = ws.atmospheric_field(50e3, 0, 0)

        pyarts.plots.AtmPoint.plot(atm_point, keys=["t", "p"])

    Parameters
    ----------
    atm_point : ~pyarts3.arts.AtmPoint
        An atmospheric point containing various atmospheric parameters
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    keys : list[str] | None, optional
        List of keys to plot. If None, plots all available keys. Defaults to None.
    **kwargs
        Additional keyword arguments passed to matplotlib plotting functions

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : list
        List of matplotlib axes objects.
    """
    if keys is None:
        keys = list(atm_point.keys())
    
    N = len(keys)
    n = int(np.ceil(np.sqrt(N)))
    
    if fig is None:
        fig = plt.figure(figsize=(4 * n, 3 * n), constrained_layout=True)
    
    if ax is None:
        ax = []
        for i, key in enumerate(keys):
            ax.append(fig.add_subplot(n, n, i + 1))
    
    for i, key in enumerate(keys):
        subplot = ax[i] if isinstance(ax, list) else ax
        value = atm_point[key]
        
        # Handle different value types
        if isinstance(value, (float, int)):
            subplot.bar([key], [value], **kwargs)
            subplot.set_ylabel("Value")
        else:
            # For array-like values, plot as line
            subplot.plot(value, **kwargs)
            subplot.set_ylabel(key)
            subplot.set_xlabel("Index")
        
        subplot.set_title(key)
        subplot.grid(True, alpha=0.3)
    
    return fig, ax
