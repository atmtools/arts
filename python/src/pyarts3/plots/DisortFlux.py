""" Plotting routine for DisortFlux """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    disort_flux: pyarts.arts.DisortFlux,
    *,
    fig=None,
    ax=None,
    freq_idx: int | None = None,
    **kwargs,
):
    """Plot DISORT flux results.

    Parameters
    ----------
    disort_flux : ~pyarts3.arts.DisortFlux
        A DisortFlux object containing flux results
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        Not used (function creates its own subplots). Accepted for API consistency.
    freq_idx : int | None, optional
        Frequency index to plot. If None, plots all frequencies. Defaults to None.
    **kwargs
        Additional keyword arguments passed to matplotlib plot/pcolormesh functions.

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : list
        List of matplotlib axes objects.
    """
    if fig is None:
        fig = plt.figure(figsize=(14, 5), constrained_layout=True)
    
    freq_grid = disort_flux.frequency_grid
    alt_grid = disort_flux.altitude_grid
    upwelling = disort_flux.up
    downwelling_diffuse = disort_flux.down_diffuse
    downwelling_direct = disort_flux.down_direct
    
    if ax is None:
        ax = []
        for i in range(3):
            ax.append(fig.add_subplot(1, 3, i + 1))
    
    if freq_idx is not None:
        # Plot single frequency
        ax[0].plot(upwelling[freq_idx, :], alt_grid, **kwargs)
        ax[0].set_xlabel('Upwelling Flux')
        ax[0].set_ylabel('Altitude [m]')
        ax[0].set_title(f'Upwelling (f={freq_grid[freq_idx]:.2e} Hz)')
        ax[0].grid(True, alpha=0.3)
        
        ax[1].plot(downwelling_diffuse[freq_idx, :], alt_grid, **kwargs)
        ax[1].set_xlabel('Diffuse Downwelling Flux')
        ax[1].set_ylabel('Altitude [m]')
        ax[1].set_title(f'Diffuse Downwelling (f={freq_grid[freq_idx]:.2e} Hz)')
        ax[1].grid(True, alpha=0.3)
        
        ax[2].plot(downwelling_direct[freq_idx, :], alt_grid, **kwargs)
        ax[2].set_xlabel('Direct Downwelling Flux')
        ax[2].set_ylabel('Altitude [m]')
        ax[2].set_title(f'Direct Downwelling (f={freq_grid[freq_idx]:.2e} Hz)')
        ax[2].grid(True, alpha=0.3)
    else:
        # Plot all frequencies as 2D
        # Use imshow for proper data/grid alignment
        # Note: alt_grid values are reversed in 'extent' ([alt_grid[-1], alt_grid[0]])
        # to ensure altitude increases upwards in the plot (matplotlib's 'imshow' by default
        # places the origin at the lower left, so this sets y-axis from lowest to highest altitude).
        extent = [freq_grid[0], freq_grid[-1], alt_grid[-1], alt_grid[0]]
        
        im1 = ax[0].imshow(upwelling.T, aspect='auto', cmap='viridis', 
                          extent=extent, origin='lower', **kwargs)
        plt.colorbar(im1, ax=ax[0], label='Flux')
        ax[0].set_xlabel('Frequency [Hz]')
        ax[0].set_ylabel('Altitude [m]')
        ax[0].set_title('Upwelling Flux')
        
        im2 = ax[1].imshow(downwelling_diffuse.T, aspect='auto', cmap='viridis',
                          extent=extent, origin='lower', **kwargs)
        plt.colorbar(im2, ax=ax[1], label='Flux')
        ax[1].set_xlabel('Frequency [Hz]')
        ax[1].set_ylabel('Altitude [m]')
        ax[1].set_title('Diffuse Downwelling Flux')
        
        im3 = ax[2].imshow(downwelling_direct.T, aspect='auto', cmap='viridis',
                          extent=extent, origin='lower', **kwargs)
        plt.colorbar(im3, ax=ax[2], label='Flux')
        ax[2].set_xlabel('Frequency [Hz]')
        ax[2].set_ylabel('Altitude [m]')
        ax[2].set_title('Direct Downwelling Flux')
    
    return fig, ax
