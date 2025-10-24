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
    alt_idx: int = 0,
    azi_idx: int = 0,
    xlabel: str = "Zenith angle [°]",
    ylabel: str = "Frequency [Hz]",
    title: str | None = None,
    **kwargs,
):
    """Plot DISORT radiance results as two 2D heatmaps (upward and downward).

    Parameters
    ----------
    disort_radiance : ~pyarts3.arts.DisortRadiance
        A DisortRadiance object containing radiance results
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        Not used (function creates its own subplots). Accepted for API consistency.
    alt_idx : int, optional
        Altitude index to use. Defaults to 0 (top of atmosphere).
    azi_idx : int, optional
        Azimuth index to use. Defaults to 0.
    xlabel : str, optional
        Label for x-axis. Defaults to "Zenith angle [°]".
    ylabel : str, optional
        Label for y-axis. Defaults to "Frequency [Hz]".
    title : str | None, optional
        Plot title. If None, generates automatic title.
    **kwargs
        Additional keyword arguments passed to matplotlib imshow function.

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : list
        List of two matplotlib axes objects [ax_up, ax_down].
    """
    if fig is None:
        fig = plt.figure(figsize=(14, 6), constrained_layout=True)
    
    freq_grid = disort_radiance.frequency_grid
    alt_grid = disort_radiance.altitude_grid
    azimuth_grid = disort_radiance.azimuth_grid
    zenith_grid = disort_radiance.zenith_grid
    radiance = disort_radiance.data
    
    # Data shape is [freq, alt, azimuth, zenith]
    # Extract slice: all frequencies, given altitude, given azimuth, all zeniths
    data_slice = radiance[:, alt_idx, azi_idx, :]
    
    # Split zenith angles: zenith > 90° is downward, zenith <= 90° is upward
    # Assuming zenith_grid is in degrees and sorted
    n_zenith = len(zenith_grid)
    mid_idx = n_zenith // 2
    
    # Upward radiation (zenith angles < 90°, smaller indices)
    upward_data = data_slice[:, :mid_idx].T  # Transpose to put freq on x-axis
    upward_zenith = zenith_grid[:mid_idx]
    
    # Downward radiation (zenith angles >= 90°, larger indices)
    downward_data = data_slice[:, mid_idx:].T  # Transpose to put freq on x-axis
    downward_zenith = zenith_grid[mid_idx:]
    
    # Create two subplots
    ax = []
    ax.append(fig.add_subplot(1, 2, 1))
    ax.append(fig.add_subplot(1, 2, 2))
    
    # Plot upward radiation
    extent_up = [freq_grid[0], freq_grid[-1], upward_zenith[0], upward_zenith[-1]]
    im1 = ax[0].imshow(upward_data, aspect='auto', cmap='viridis', 
                       extent=extent_up, origin='lower', **kwargs)
    plt.colorbar(im1, ax=ax[0], label='Radiance')
    ax[0].set_xlabel(ylabel)  # Swapped
    ax[0].set_ylabel(xlabel)  # Swapped
    ax[0].set_title('Upward Radiation')
    
    # Plot downward radiation
    extent_down = [freq_grid[0], freq_grid[-1], downward_zenith[0], downward_zenith[-1]]
    im2 = ax[1].imshow(downward_data, aspect='auto', cmap='viridis', 
                       extent=extent_down, origin='lower', **kwargs)
    plt.colorbar(im2, ax=ax[1], label='Radiance')
    ax[1].set_xlabel(ylabel)  # Swapped
    ax[1].set_ylabel(xlabel)  # Swapped
    ax[1].set_title('Downward Radiation')
    
    if title is not None:
        fig.suptitle(title)
    else:
        fig.suptitle(f'DISORT Radiance (alt={alt_grid[alt_idx]:.1f} m, φ={azimuth_grid[azi_idx]:.1f}°)')
    
    return fig, ax
