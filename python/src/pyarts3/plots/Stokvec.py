""" Plotting routine for Stokvec """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    stokvec: pyarts.arts.Stokvec,
    *,
    fig=None,
    ax=None,
    labels: list[str] = ['I', 'Q', 'U', 'V'],
    title: str = "Stokes Vector",
):
    """Plot a Stokes vector as a bar chart.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts

        # Create a Stokes vector (unpolarized radiation)
        stokes = pyarts.arts.Stokvec([1.0, 0.1, 0.05, 0.0])

        pyarts.plots.Stokvec.plot(stokes)

    Parameters
    ----------
    stokvec : ~pyarts3.arts.Stokvec
        A Stokes vector (4 components: I, Q, U, V)
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    labels : list[str], optional
        Labels for the 4 Stokes parameters. Defaults to ['I', 'Q', 'U', 'V'].
    title : str, optional
        Plot title. Defaults to "Stokes Vector".

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : As input
        The matplotlib axes.
    """
    if fig is None:
        fig = plt.figure(figsize=(8, 6))
    
    if ax is None:
        ax = fig.add_subplot(1, 1, 1)
    
    ax.bar(labels, stokvec)
    ax.set_ylabel("Value")
    ax.set_title(title)
    ax.grid(True, alpha=0.3, axis='y')
    
    return fig, ax
