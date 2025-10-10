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
    freqs: np.ndarray | float = None,
    is_polar: bool | None = None,
):
    """Plot the transmission matrix elements.

    Selects a 1D plot if the transmission matrix is 1D, otherwise a 4x4 grid of plots.

    Parameters
    ----------
    transmission_matrix : ~pyarts3.arts.MuelmatVector
        A vector of Mueller matrices representing the transmission properties.
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
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
    subs : As input
        As input.
    """
    if is_polar is None:
        is_polar = transmission_matrix.is_polarized()
    
    if freqs is None:
        freqs = np.arange(len(transmission_matrix))
    if len(freqs) != len(transmission_matrix):
        raise ValueError("Length of freqs must match length of transmission_matrix")
    
    if fig is None:
        if not is_polar:
            fig, subs = plt.subplots(1, 1, figsize=(6, 4))
        else:
            fig, subs = plt.subplots(4, 4, figsize=(12, 12), constrained_layout=True)
    else:
        subs = fig.axes
        if not is_polar and len(subs) != 1:
            raise ValueError("Figure must have exactly one axis for 1D transmission matrix")
        elif len(subs) != 16:
            raise ValueError("Figure must have exactly 16 axes for 4x4 transmission matrix")

    if not is_polar:
        subs.plot(freqs, [m[0, 0] for m in transmission_matrix], label='I')
    else:
        for i in range(4):
            for j in range(4):
                ax = subs[i, j]
                ax.plot(freqs, [m[i, j] for m in transmission_matrix])

    return fig, subs
