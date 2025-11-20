""" Plotting routine for DisortRadiance """

import numpy
import matplotlib
import pyarts3 as pyarts
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.DisortRadiance,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         alt_idx: int = 0,
         azi_idx: int = 0,
         plotstyle='contourf',
         select: list = ['up', 'down'],
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Plot DISORT radiance results as two 2D heatmaps (upward and downward).

    .. rubric:: Example (mock data)

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        radiance = pyarts.arts.DisortRadiance()
        radiance.freq_grid = pyarts.arts.AscendingGrid(np.linspace(1e9, 1e12, 50))
        radiance.alt_grid = pyarts.arts.DescendingGrid(np.linspace(80e3, 0, 41))
        radiance.azi_grid = pyarts.arts.AziGrid(np.linspace(0, 330, 12))
        radiance.zen_grid = pyarts.arts.ZenGrid(np.linspace(0, 180, 17))
        radiance.data = pyarts.arts.Tensor4(np.random.randn(50, 41, 12, 17) * 10 + 100 * np.sin(np.radians(radiance.zen_grid))[None, None, None, :])
        fig, ax = pyarts.plots.DisortRadiance.plot(radiance)

    Parameters
    ----------
    data : ~pyarts3.arts.DisortRadiance
        A DisortRadiance object containing radiance results
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes] | None, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    alt_idx : int, optional
        Altitude index to use. Defaults to 0 (top of atmosphere).
    azi_idx : int, optional
        Azimuth index to use. Defaults to 0.
    plotstyle : str, optional
        The matplotlib plotting function to use, 'contourf' or 'plot'.
    select : list or str, optional
        Which directions to plot: 'up' and/or 'down'. Defaults to ['up', 'down'].
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig :
        As input if input.  Otherwise the created Figure.
    ax :
        As input if input.  Otherwise the created Axes.
    """
    freq_grid = data.freq_grid
    zen_grid = data.zen_grid
    radiance = data.data

    # Data shape is [freq, alt, azimuth, zenith]
    # Extract slice: all frequencies, given altitude, given azimuth, all zeniths
    data_slice = radiance[:, alt_idx, azi_idx, :]

    # Split zenith angles: zenith > 90° is downward, zenith <= 90° is upward
    # Assuming zen_grid is in degrees and sorted
    n_zenith = len(zen_grid)
    mid_idx = n_zenith // 2

    # Upward radiation (zenith angles < 90°, smaller indices)
    upward_data = data_slice[:, mid_idx:]  # Transpose to put freq on x-axis
    upward_zenith = zen_grid[mid_idx:]

    # Downward radiation (zenith angles >= 90°, larger indices)
    downward_data = data_slice[:, :mid_idx]  # Transpose to put freq on x-axis
    downward_zenith = zen_grid[:mid_idx]

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
            select_flat_ax(ax, 0).contourf(
                freq_grid, upward_zenith, upward_data.T, **kwargs)

        # Plot downward radiation
        if has_down:
            select_flat_ax(ax, n - 1).contourf(freq_grid,
                                               downward_zenith, downward_data.T, **kwargs)
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
