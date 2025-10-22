""" Plotting routine for DisortRadiance """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    disort_radiance: pyarts.arts.DisortRadiance,
    *,
    fig=None,
    ax=None,
    freq_idx: int = 0,
    alt_idx: int = 0,
    xlabel: str = "Azimuth [°]",
    ylabel: str = "Zenith [°]",
    title: str | None = None,
    colorbar: bool = True,
    cmap: str = "viridis",
    **kwargs,
):
    """Plot DISORT radiance results.

    Parameters
    ----------
    disort_radiance : ~pyarts3.arts.DisortRadiance
        A DisortRadiance object containing radiance results
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    freq_idx : int, optional
        Frequency index to plot. Defaults to 0.
    alt_idx : int, optional
        Altitude index to plot. Defaults to 0.
    xlabel : str, optional
        Label for x-axis. Defaults to "Azimuth [°]".
    ylabel : str, optional
        Label for y-axis. Defaults to "Zenith [°]".
    title : str | None, optional
        Plot title. If None, generates automatic title. Defaults to None.
    colorbar : bool, optional
        Whether to show colorbar. Defaults to True.
    cmap : str, optional
        Colormap name. Defaults to "viridis".
    **kwargs
        Additional keyword arguments passed to matplotlib pcolormesh function.

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : As input
        The matplotlib axes.
    """
    if fig is None:
        fig = plt.figure(figsize=(10, 8))
    
    if ax is None:
        ax = fig.add_subplot(1, 1, 1)
    
    freq_grid = disort_radiance.frequency_grid
    alt_grid = disort_radiance.level_altitude_grid
    azimuth_grid = disort_radiance.azimuth_grid
    zenith_grid = disort_radiance.zenith_grid
    radiance = disort_radiance.radiance_data
    
    # Extract slice for given frequency and altitude
    data_slice = radiance[freq_idx, alt_idx, :, :]
    
    # Create meshgrid for plotting
    zen_mesh, azi_mesh = np.meshgrid(zenith_grid, azimuth_grid, indexing='ij')
    
    im = ax.pcolormesh(azi_mesh, zen_mesh, data_slice, cmap=cmap, shading='auto', **kwargs)
    
    if colorbar:
        plt.colorbar(im, ax=ax, label='Radiance')
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    if title is None:
        title = f'DISORT Radiance (f={freq_grid[freq_idx]:.2e} Hz, alt={alt_grid[alt_idx]:.1f} m)'
    ax.set_title(title)
    
    return fig, ax
