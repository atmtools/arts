""" Plotting routine for Sun """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    sun: pyarts.arts.Sun,
    *,
    fig=None,
    ax=None,
    freqs: np.ndarray | None = None,
    xlabel: str = "Frequency [Hz]",
    ylabel: str = "Spectrum",
    title: str = "Solar Spectrum",
    show_info: bool = True,
    **kwargs
):
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
    sun : ~pyarts3.arts.Sun
        A Sun object containing spectrum and properties
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    freqs : :class:`~numpy.ndarray` | None, optional
        Frequency grid for x-axis. If None, uses indices. Defaults to None.
    xlabel : str, optional
        Label for x-axis. Defaults to "Frequency [Hz]".
    ylabel : str, optional
        Label for y-axis. Defaults to "Spectrum".
    title : str, optional
        Plot title. Defaults to "Solar Spectrum".
    show_info : bool, optional
        Whether to display info text box with sun properties. Defaults to True.
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
    
    # spectrum is a Matrix (nfreq x 4 Stokes), use first Stokes dimension
    spectrum = sun.spectrum[:, 0] if sun.spectrum.shape[1] > 0 else sun.spectrum[:, 0]
    
    if freqs is None:
        freqs = np.arange(len(spectrum))
    
    ax.plot(freqs, spectrum, **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    
    # Add text with sun properties
    if show_info:
        info_text = f"Radius: {sun.radius:.2e} m\n"
        info_text += f"Distance: {sun.distance:.2e} m\n"
        info_text += f"Latitude: {sun.latitude:.2f}°\n"
        info_text += f"Longitude: {sun.longitude:.2f}°"
        
        ax.text(0.02, 0.98, info_text, transform=ax.transAxes,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    return fig, ax
