""" Plotting routine for SortedGriddedField1 """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(
    gridded_field: pyarts.arts.SortedGriddedField1,
    *,
    fig=None,
    ax=None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str | None = None,
    **kwargs
):
    """Plot a SortedGriddedField1 as a line plot using its grid.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        # Create a gridded field
        freqs = pyarts.arts.AscendingGrid(np.linspace(1e9, 1e12, 100))
        data = pyarts.arts.Vector(np.exp(-((freqs - 5e11)/1e11)**2))
        
        field = pyarts.arts.SortedGriddedField1()
        field.grids = (freqs,)
        field.data = data
        field.gridnames = ("Frequency [Hz]",)
        field.dataname = "Response"

        pyarts.plots.SortedGriddedField1.plot(field)

    Parameters
    ----------
    gridded_field : ~pyarts3.arts.SortedGriddedField1
        A 1D sorted gridded field with named grid
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    xlabel : str | None, optional
        Label for x-axis. If None, uses grid name. Defaults to None.
    ylabel : str | None, optional
        Label for y-axis. If None, uses data name. Defaults to None.
    title : str | None, optional
        Plot title. If None, uses data name. Defaults to None.
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
        fig = plt.figure(figsize=(8, 6))
    
    if ax is None:
        ax = fig.add_subplot(1, 1, 1)
    
    # Get grid and data
    xgrid = gridded_field.grids[0]
    data = gridded_field.data
    
    ax.plot(xgrid, data, **kwargs)
    ax.set_xlabel(xlabel if xlabel is not None else gridded_field.gridnames[0])
    ax.set_ylabel(ylabel if ylabel is not None else gridded_field.dataname)
    ax.set_title(title if title is not None else gridded_field.dataname)
    ax.grid(True, alpha=0.3)
    
    return fig, ax
