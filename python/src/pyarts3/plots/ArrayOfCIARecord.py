from pyarts3.arts import ArrayOfCIARecord, AscendingGrid, AtmPoint, AtmField, PropagationPathPoint
import pyarts3 as pyarts
from .common import default_fig_ax, select_flat_ax
from . import CIARecord

__all__ = [
    'plot',
]


def plot(data: ArrayOfCIARecord,
         *,
         fig=None,
         ax=None,
         same: bool = False,
         freqs: AscendingGrid = None,
         atm: AtmPoint | AtmField | None = None,
         path_point: PropagationPathPoint = PropagationPathPoint(),
         T_extrapolfac: float = 1e99,
         **kwargs):
    """Creates a plot of the absorption by the records in the ArrayOfCIARecord object.

    The plot passes all relevant arguments to the CIARecord plotting routine.

    If freqs is given, this range will be used, otherwise the built-in frequency range
    is used.

    If atm is None, an atmospheric point is created by reading the tropical standard atmosphere
    from the ARTS database and extracting the AtmPoint at position by the path_point.  If atm is
    an AtmField, the AtmPoint at position by the path_point is extracted.

    The path_point is passed directly to the ArrayOfCIARecord.spectral_propmat() method and defaults
    to pos: [0, 0, 0], los: [0, 0].

    same determines if each entry in the array is plotted onto the same plot or not.

    Parameters
    ----------
    data : ~pyarts3.arts.ArrayOfCIARecord
        The ArrayOfCIARecord object containing the data to plot.
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    same : bool, optional
        Draw on a single canvas or not.
    freqs : ~pyarts3.arts.AscendingGrid | int, optional
        The frequency grid to use for the x-axis. Defaults to None manually creating a range.
    atm : ~pyarts3.arts.AtmPoint | ~pyarts3.arts.AtmField | None, optional
        The atmospheric point data to use for the plot. Defaults to None computed values.
    path_point : ~pyarts3.arts.PropagationPathPoint, optional
        The propagation path point to use for the plot. Defaults to pos: [0, 0, 0], los: [0, 0].
    T_extrapolfac : float, optional
        Internal extrapolation factor
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.  Ignores label.

    Returns
    -------
    fig : As input
        As input.
    ax : As input
        As input.
    """

    if atm is None:
        ws = pyarts.Workspace()
        ws.abs_speciesSet(species=[f"{band.isot}" for band in data])
        basename = "planets/Earth/afgl/tropical/"
        toa = 1 + path_point.pos[0]
        ws.atm_fieldRead(toa=toa, basename=basename, missing_is_zero=1)
        atm = ws.atm_field(*path_point.pos)
    elif isinstance(atm, AtmField):
        atm = atm(*path_point.pos)

    n = 1 if same else len(data)
    fig, ax = default_fig_ax(fig, ax, n, 1,  fig_kwargs={'figsize': (6, 4*n)})

    if "label" in kwargs:
        del kwargs["label"]

    for i, v in enumerate(data):
        CIARecord.plot(v, fig=fig, ax=select_flat_ax(ax, 0 if same else i),
                        freqs=freqs, atm=atm, path_point=path_point,
                        T_extrapolfac=T_extrapolfac, label=f"{v.specs[0]}-CIA-{v.specs[1]}",
                        **kwargs)

    return fig, ax
