""" Plotting routine for DisortRadiance """

import pyarts3 as pyarts
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(
    data: pyarts.arts.DisortRadiance,
    *,
    fig=None,
    ax=None,
    alt_idx: int = 0,
    azi_idx: int = 0,
    plotstyle='contourf',
    select: list = ['up', 'down'],
    freqs=None,
    **kwargs,
):
    """Plot DISORT radiance results as two 2D heatmaps (upward and downward).

    Parameters
    ----------
    data : ~pyarts3.arts.DisortRadiance
        A DisortRadiance object containing radiance results
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        Not used (function creates its own subplots). Accepted for API consistency.
    alt_idx : int, optional
        Altitude index to use. Defaults to 0 (top of atmosphere).
    azi_idx : int, optional
        Azimuth index to use. Defaults to 0.
    plotstyle : str, optional
        The matplotlib plotting function to use, 'contourf' or 'plot'
    freqs : array-like, optional
        Frequency grid to use for x-axis. Defaults to None, which uses data.frequency_grid.
    select : list or str, optional
        Which directions to plot: 'up' and/or 'down'. Defaults to ['up', 'down'].
    **kwargs
        Additional keyword arguments passed to matplotlib plotting function.

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : list
        List of two matplotlib axes objects [ax_up, ax_down].
    """
    freq_grid = data.frequency_grid if freqs is None else freqs
    zenith_grid = data.zenith_grid
    radiance = data.data

    # Data shape is [freq, alt, azimuth, zenith]
    # Extract slice: all frequencies, given altitude, given azimuth, all zeniths
    data_slice = radiance[:, alt_idx, azi_idx, :]

    # Split zenith angles: zenith > 90° is downward, zenith <= 90° is upward
    # Assuming zenith_grid is in degrees and sorted
    n_zenith = len(zenith_grid)
    mid_idx = n_zenith // 2

    # Upward radiation (zenith angles < 90°, smaller indices)
    upward_data = data_slice[:, mid_idx:]  # Transpose to put freq on x-axis
    upward_zenith = zenith_grid[mid_idx:]

    # Downward radiation (zenith angles >= 90°, larger indices)
    downward_data = data_slice[:, :mid_idx]  # Transpose to put freq on x-axis
    downward_zenith = zenith_grid[:mid_idx]

    # Plotting options
    select = [select] if isinstance(select, str) else select
    has_up = 'up' in select
    has_down = 'down' in select
    n = has_up + has_down
    if n == 0:
        raise ValueError(f"select option unknown '{select}'. See docs.")

    fig, ax = default_fig_ax(fig, ax, 1, n, fig_kwargs={
                             'figsize': (n*7, 6), 'constrained_layout': True})

    if plotstyle == 'contourf':
        # Plot upward radiation
        if has_up:
            select_flat_ax(ax, 0).contourf(freq_grid, upward_zenith, upward_data.T, **kwargs)

        # Plot downward radiation
        if has_down:
            select_flat_ax(ax, n - 1).contourf(freq_grid, downward_zenith, downward_data.T, **kwargs)
    elif plotstyle == 'plot':
        # Plot upward radiation
        if has_up:
            select_flat_ax(ax, 0).plot(freq_grid, upward_data, **kwargs)

        # Plot downward radiation
        if has_down:
            select_flat_ax(ax, n - 1).plot(freq_grid, downward_data, **kwargs)
    else:
        raise ValueError(f"plotstyle unknown '{plotstyle}'. See docs.")

    return fig, ax
