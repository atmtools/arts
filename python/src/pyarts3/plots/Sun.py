""" Plotting routine for Sun """

import numpy
import matplotlib
import pyarts3 as pyarts
import numpy as np
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.Sun,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         freqs: pyarts.arts.AscendingGrid | None = None,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Plot the solar spectrum.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create a simplified sun object
        sun = pyarts.arts.Sun()
        sun.spectrum = pyarts.arts.Vector(np.ones(100) * 1e-3)  # Accepts Vector
        sun.radius = 6.96e8  # meters
        sun.distance = 1.496e11  # meters (1 AU)
        sun.latitude = 0.0
        sun.longitude = 0.0

        pyarts.plots.Sun.plot(sun)

    Parameters
    ----------
    data : ~pyarts3.arts.Sun
        A Sun object containing spectrum and properties
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes], optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    freqs : ~pyarts3.arts.AscendingGrid | None, optional
        Frequency grid for x-axis. If None, uses indices. Defaults to None.
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig :
        As input if input.  Otherwise the created Figure.
    ax :
        As input if input.  Otherwise the created Axes.
    """
    fig, ax = default_fig_ax(fig, ax, fig_kwargs={'figsize': (10, 6)})

    # spectrum is a Matrix (nfreq x 4 Stokes), use first Stokes dimension
    spectrum = data.spectrum[:, 0]

    freqs = np.arange(len(spectrum)) if freqs is None else freqs

    select_flat_ax(ax, 0).plot(freqs, spectrum, **kwargs)

    return fig, ax
