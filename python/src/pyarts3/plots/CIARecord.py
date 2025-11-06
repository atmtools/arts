from pyarts3.arts import CIARecord, AscendingGrid, AtmPoint, AtmField, Propmat, SpeciesEnum, PropagationPathPoint
import pyarts3 as pyarts
import matplotlib
import numpy as np
from .common import default_fig_ax, select_flat_ax
from copy import deepcopy as copy

__all__ = [
    'plot',
]


def plot(data: CIARecord,
         *,
         fig=None,
         ax=None,
         freqs: AscendingGrid = None,
         atm: AtmPoint | AtmField | None = None,
         path_point: PropagationPathPoint = PropagationPathPoint(),
         T_extrapolfac: float = 1e99,
         **kwargs):
    """Creates a plot of the absorption by the records in the CIARecord object.

    If freqs is given, this range will be used, otherwise the built-in frequency range
    is used.

    If atm is None, and atmospheric point is created by reading the tropical standard atmosphere
    from the ARTS database and extracting the AtmPoint at position by the path_point.  If atm is
    an AtmField, the AtmPoint at position by the path_point is extracted.

    The path_point is passed directly to the CIARecord.propagation_matrix() method and defaults
    to pos: [0, 0, 0], los: [0, 0].

    Parameters
    ----------
    data : ~pyarts3.arts.CIARecord
        The CIARecord object containing the data to plot.
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    freqs : ~pyarts3.arts.AscendingGrid | int, optional
        The frequency grid to use for the x-axis. Defaults to None manually creating a range.
    atm : ~pyarts3.arts.AtmPoint | ~pyarts3.arts.AtmField | None, optional
        The atmospheric point data to use for the plot. Defaults to None computed values.
    path_point : ~pyarts3.arts.PropagationPathPoint, optional
        The propagation path point to use for the plot. Defaults to pos: [0, 0, 0], los: [0, 0].
    T_extrapolfac : float, optional
        Internal extrapolation factor
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig : As input
        As input.
    ax : As input
        As input.
    """

    if atm is None:
        ws = pyarts.Workspace()
        ws.absorption_speciesSet(species=[f"{band.isot}" for band in data])
        basename = "planets/Earth/afgl/tropical/"
        toa = 1 + path_point.pos[0]
        ws.atmospheric_fieldRead(toa=toa, basename=basename, missing_is_zero=1)
        atm = ws.atmospheric_field(*path_point.pos)
    elif isinstance(atm, AtmField):
        atm = atm(*path_point.pos)

    fig, ax = default_fig_ax(fig, ax, 1,  1, fig_kwargs={'figsize': (6, 4)})

    if freqs is None:
        for i in range(len(data.data)):
            cia = copy(data)
            cia.data = [data.data[i]]
            freqs = cia.data[0].grids[0]
            select_flat_ax(ax, 0).plot(
                freqs,
                cia.propagation_matrix(f=freqs, atm=atm, T_extrapolfac=T_extrapolfac)[:, 0],
                **kwargs
            )
    else:
        select_flat_ax(ax, 0).plot(
            freqs,
            data.propagation_matrix(f=freqs, atm=atm, T_extrapolfac=T_extrapolfac)[:, 0],
            **kwargs
        )

    return fig, ax
