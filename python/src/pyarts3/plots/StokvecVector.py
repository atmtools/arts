""" Plotting routine for StokvecVector """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    stokvec_vector: pyarts.arts.StokvecVector,
    *,
    fig=None,
    ax=None,
    freqs: np.ndarray | None = None,
    xlabel: str = "Index",
    ylabel: str = "Value",
    title: str = "Stokes Vector Components",
    labels: list[str] = ['I', 'Q', 'U', 'V'],
    **kwargs
):
    """Plot a vector of Stokes vectors showing all 4 components.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create spectral radiance (Stokes vectors at each frequency)
        n = 100
        stokes_vec = pyarts.arts.StokvecVector(np.zeros((n, 4)))
        stokes_vec[:, 0] = 1 + 0.5 * np.sin(np.linspace(0, 4*np.pi, n))  # I
        stokes_vec[:, 1] = 0.1 * np.cos(np.linspace(0, 4*np.pi, n))      # Q

        pyarts.plots.StokvecVector.plot(stokes_vec, xlabel="Frequency Index")

    Parameters
    ----------
    stokvec_vector : ~pyarts3.arts.StokvecVector
        A vector of Stokes vectors (each with 4 components: I, Q, U, V)
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        Not used (function creates its own subplots). Accepted for API consistency.
    freqs : :class:`~numpy.ndarray` | None, optional
        Frequency or position grid for x-axis. If None, uses indices. Defaults to None.
    xlabel : str, optional
        Label for x-axis. Defaults to "Index".
    ylabel : str, optional
        Label for y-axis. Defaults to "Value".
    title : str, optional
        Plot title. Defaults to "Stokes Vector Components".
    labels : list[str], optional
        Labels for the 4 Stokes parameters. Defaults to ['I', 'Q', 'U', 'V'].
    **kwargs
        Additional keyword arguments passed to matplotlib plot()

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : list
        List of matplotlib axes objects (4 subplots).
    """
    if fig is None:
        fig = plt.figure(figsize=(12, 10), constrained_layout=True)
    
    if ax is None:
        ax = []
        for i in range(4):
            ax.append(fig.add_subplot(2, 2, i + 1))
    
    if freqs is None:
        freqs = np.arange(len(stokvec_vector))
    
    for i, label in enumerate(labels):
        subplot = ax[i] if isinstance(ax, list) else ax
        data = stokvec_vector[:, i]
        subplot.plot(freqs, data, label=label, **kwargs)
        subplot.set_xlabel(xlabel)
        subplot.set_ylabel(ylabel)
        subplot.set_title(f"{title} - {label}")
        subplot.legend()
        subplot.grid(True, alpha=0.3)
    
    return fig, ax
