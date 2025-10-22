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
    freq_idx: int | None = None,
):
    """Plot DISORT flux results.

    Parameters
    ----------
    disort_flux : ~pyarts3.arts.DisortFlux
        A DisortFlux object containing flux results
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    freq_idx : int | None, optional
        Frequency index to plot. If None, plots all frequencies. Defaults to None.

    Returns
    -------
    fig : As input
        The matplotlib figure.
    axes : list
        List of matplotlib axes objects.
    """
    if fig is None:
        fig = plt.figure(figsize=(14, 5), constrained_layout=True)
    
    freq_grid = disort_flux.frequency_grid
    alt_grid = disort_flux.level_altitude_grid
    upwelling = disort_flux.upwelling_flux
    downwelling_diffuse = disort_flux.diffuse_downwelling_flux
    downwelling_direct = disort_flux.direct_downwelling_flux
    
    axes = []
    
    if freq_idx is not None:
        # Plot single frequency
        ax1 = fig.add_subplot(1, 3, 1)
        ax1.plot(upwelling[freq_idx, :], alt_grid)
        ax1.set_xlabel('Upwelling Flux')
        ax1.set_ylabel('Altitude [m]')
        ax1.set_title(f'Upwelling (f={freq_grid[freq_idx]:.2e} Hz)')
        ax1.grid(True, alpha=0.3)
        axes.append(ax1)
        
        ax2 = fig.add_subplot(1, 3, 2)
        ax2.plot(downwelling_diffuse[freq_idx, :], alt_grid)
        ax2.set_xlabel('Diffuse Downwelling Flux')
        ax2.set_ylabel('Altitude [m]')
        ax2.set_title(f'Diffuse Downwelling (f={freq_grid[freq_idx]:.2e} Hz)')
        ax2.grid(True, alpha=0.3)
        axes.append(ax2)
        
        ax3 = fig.add_subplot(1, 3, 3)
        ax3.plot(downwelling_direct[freq_idx, :], alt_grid)
        ax3.set_xlabel('Direct Downwelling Flux')
        ax3.set_ylabel('Altitude [m]')
        ax3.set_title(f'Direct Downwelling (f={freq_grid[freq_idx]:.2e} Hz)')
        ax3.grid(True, alpha=0.3)
        axes.append(ax3)
    else:
        # Plot all frequencies as 2D
        ax1 = fig.add_subplot(1, 3, 1)
        im1 = ax1.pcolormesh(freq_grid, alt_grid, upwelling.T, cmap='viridis', shading='auto')
        plt.colorbar(im1, ax=ax1, label='Flux')
        ax1.set_xlabel('Frequency [Hz]')
        ax1.set_ylabel('Altitude [m]')
        ax1.set_title('Upwelling Flux')
        axes.append(ax1)
        
        ax2 = fig.add_subplot(1, 3, 2)
        im2 = ax2.pcolormesh(freq_grid, alt_grid, downwelling_diffuse.T, cmap='viridis', shading='auto')
        plt.colorbar(im2, ax=ax2, label='Flux')
        ax2.set_xlabel('Frequency [Hz]')
        ax2.set_ylabel('Altitude [m]')
        ax2.set_title('Diffuse Downwelling Flux')
        axes.append(ax2)
        
        ax3 = fig.add_subplot(1, 3, 3)
        im3 = ax3.pcolormesh(freq_grid, alt_grid, downwelling_direct.T, cmap='viridis', shading='auto')
        plt.colorbar(im3, ax=ax3, label='Flux')
        ax3.set_xlabel('Frequency [Hz]')
        ax3.set_ylabel('Altitude [m]')
        ax3.set_title('Direct Downwelling Flux')
        axes.append(ax3)
    
    return fig, axes
