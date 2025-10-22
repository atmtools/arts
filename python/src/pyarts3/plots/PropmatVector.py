""" Plotting routine for PropmatVector """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    propmat_vector: pyarts.arts.PropmatVector,
    *,
    fig=None,
    ax=None,
    freqs: np.ndarray | None = None,
    component: int = 0,
    xlabel: str = "Index",
    ylabel: str = "Propagation Matrix Component",
    title: str = "Propagation Matrix Vector",
    **kwargs
):
    """Plot a propagation matrix vector.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create a propagation matrix vector (absorption spectrum)
        n = 50
        propmat_vec = pyarts.arts.PropmatVector(np.zeros((n, 7)))
        propmat_vec[:, 0] = np.exp(-((np.arange(n) - 25)/10)**2)  # Absorption

        pyarts.plots.PropmatVector.plot(propmat_vec, component=0,
                                        ylabel="Absorption [1/m]")

    Parameters
    ----------
    propmat_vector : ~pyarts3.arts.PropmatVector
        A vector of propagation matrices (7 components each)
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    freqs : :class:`~numpy.ndarray` | None, optional
        Frequency or position grid for x-axis. If None, uses indices. Defaults to None.
    component : int, optional
        Which component to plot (0-6). Defaults to 0.
    xlabel : str, optional
        Label for x-axis. Defaults to "Index".
    ylabel : str, optional
        Label for y-axis. Defaults to "Propagation Matrix Component".
    title : str, optional
        Plot title. Defaults to "Propagation Matrix Vector".
    **kwargs
        Additional keyword arguments passed to plot()

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : As input
        The matplotlib axes.
    """
    if fig is None:
        fig = plt.figure(figsize=(10, 6))
    
    if ax is None:
        ax = fig.add_subplot(1, 1, 1)
    
    # Extract component from propmat vector
    data = propmat_vector[:, component]
    
    if freqs is None:
        freqs = np.arange(len(data))
    
    ax.plot(freqs, data, **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"{title} (Component {component})")
    ax.grid(True, alpha=0.3)
    
    return fig, ax
