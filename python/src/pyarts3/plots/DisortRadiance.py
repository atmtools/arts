""" Plotting routine for DisortRadiance """

import pyarts3 as pyarts
from .common import default_fig_ax, select_flat_ax

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
    **kwargs
        Additional keyword arguments passed to matplotlib imshow function.

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : list
        List of two matplotlib axes objects [ax_up, ax_down].
    """
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

    fig, ax = default_fig_ax(fig, ax, 1, 2, fig_kwargs={
                             'figsize': (14, 6), 'constrained_layout': True})

    # Plot upward radiation
    extent_up = [freq_grid[0], freq_grid[-1], upward_zenith[0], upward_zenith[-1]]
    select_flat_ax(ax, 0).imshow(upward_data,
                                 extent=extent_up, **kwargs)

    # Plot downward radiation
    extent_down = [freq_grid[0], freq_grid[-1], downward_zenith[0], downward_zenith[-1]]
    select_flat_ax(ax, 1).imshow(downward_data,
                                 extent=extent_down, **kwargs)

    return fig, ax
