""" Plotting routine for SortedGriddedField2 """

import pyarts3 as pyarts
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(
    gridded_field: pyarts.arts.SortedGriddedField2,
    *,
    fig=None,
    ax=None,
    **kwargs
):
    """Plot a SortedGriddedField2 as a 2D heatmap using its sorted grids.

    Parameters
    ----------
    gridded_field : ~pyarts3.arts.SortedGriddedField2
        A 2D sorted gridded field with ascending grids
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    **kwargs
        Additional keyword arguments passed to pcolormesh()

    Returns
    -------
    fig : As input
        The matplotlib figure.
    ax : As input
        The matplotlib axes.
    """
    fig, ax = default_fig_ax(fig, ax, fig_kwargs={'figsize': (10, 6)})

    # Get grids and data
    xgrid = gridded_field.grids[0]
    ygrid = gridded_field.grids[1]
    data = gridded_field.data

    # Plot using pcolormesh
    select_flat_ax(ax, 0).pcolormesh(xgrid, ygrid, data, **kwargs)

    return fig, ax
