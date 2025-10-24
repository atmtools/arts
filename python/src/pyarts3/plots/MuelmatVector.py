""" Plotting routine for transmission matrix """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    transmission_matrix: pyarts.arts.MuelmatVector,
    *,
    fig=None,
    ax=None,
    freqs: np.ndarray | float = None,
    is_polar: bool | None = None,
    **kwargs,
):
    """Plot the transmission matrix elements.

    Selects a 1D plot if the transmission matrix is 1D, otherwise a 4x4 grid of plots.

    Parameters
    ----------
    transmission_matrix : ~pyarts3.arts.MuelmatVector
        A vector of Mueller matrices representing the transmission properties.
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes or list, optional
        The matplotlib axes to draw on. For 1D plot: single axes. For 4x4 plot: 16 axes.
        Defaults to None for new axes.
    freqs : list, optional
        A list of frequencies to plot. Defaults to None for no frequency grid.
    is_polar : bool, optional
        If True, the transmission matrix is treated as polarized (4x4).
        If False, it is treated as unpolarized (1D).
        If None (default), the polarization state is determined from the data.

    Returns
    -------
    fig : As input
        As input.
    ax : As input
        As input. Single axis for 1D plot, or list of 16 axes for 4x4 plot.
    """
    if is_polar is None:
        is_polar = transmission_matrix.is_polarized()
    
    if freqs is None:
        freqs = np.arange(len(transmission_matrix))
    if len(freqs) != len(transmission_matrix):
        raise ValueError("Length of freqs must match length of transmission_matrix")
    
    if fig is None:
        if not is_polar:
            fig = plt.figure(figsize=(6, 4))
        else:
            fig = plt.figure(figsize=(12, 12), constrained_layout=True)
    
    if ax is None:
        if not is_polar:
            ax = fig.add_subplot(1, 1, 1)
        else:
            ax = []
            for i in range(4):
                row = []
                for j in range(4):
                    row.append(fig.add_subplot(4, 4, i*4 + j + 1))
                ax.append(row)

    if not is_polar:
        ax.plot(freqs, [m[0, 0] for m in transmission_matrix], label='I', **kwargs)
    else:
        for i in range(4):
            for j in range(4):
                ax[i][j].plot(freqs, [m[i, j] for m in transmission_matrix], **kwargs)

    return fig, ax
