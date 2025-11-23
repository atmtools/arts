import numpy
import matplotlib
from pyarts3.arts import CIARecord, AscendingGrid, AtmPoint, AtmField, PropagationPathPoint, SpeciesEnum
import pyarts3 as pyarts
from .common import default_fig_ax, select_flat_ax
from copy import deepcopy as copy

__all__ = [
    'plot',
]


def plot(data: CIARecord,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         freqs: AscendingGrid = None,
         atm: AtmPoint | AtmField | None = None,
         spec1: SpeciesEnum = SpeciesEnum("AIR"),
         spec2: SpeciesEnum = SpeciesEnum("AIR"),
         path_point: PropagationPathPoint = PropagationPathPoint(),
         T_extrapolfac: float = 1e99,
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Creates a plot of the absorption by the records in the CIARecord object.

    If freqs is given, this range will be used, otherwise the built-in frequency range
    is used.

    If atm is None, an atmospheric point is created by reading the tropical standard atmosphere
    from the ARTS database and extracting the AtmPoint at position by the path_point.  If atm is
    an AtmField, the AtmPoint at position by the path_point is extracted.

    The path_point is passed directly to the CIARecord.spectral_propmat() method and defaults
    to pos: [0, 0, 0], los: [0, 0].

    .. rubric:: Example

    .. plot::
        :include-source:

        import matplotlib.pyplot as plt
        import pyarts3 as pyarts

        cia = pyarts.arts.CIARecord.fromxml("cia/O2-CIA-N2.xml")
        f, a = pyarts.plots.CIARecord.plot(cia, fig=plt.figure(figsize=(12, 6)), spec1="O2", spec2="N2")
        a.set_yscale("log")
        a.set_xlabel("Frequency [Hz]")
        a.set_ylabel("Absorption [1/m]")
        a.set_title("O$_2$-N$_2$ Collision-induced absorption")

    Parameters
    ----------
    data : ~pyarts3.arts.CIARecord
        The CIARecord object containing the data to plot.
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes] | None, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    freqs : ~pyarts3.arts.AscendingGrid | int, optional
        The frequency grid to use for the x-axis. Defaults to None manually creating a range.
    atm : ~pyarts3.arts.AtmPoint | ~pyarts3.arts.AtmField | None, optional
        The atmospheric point data to use for the plot. Defaults to None computed values.
    spec1 : ~pyarts3.arts.SpeciesEnum, optional
        The first species for the CIA calculation.  Defaults to "AIR" for no species.
    spec2 : ~pyarts3.arts.SpeciesEnum, optional
        The second species for the CIA calculation.  Defaults to "AIR" for no species.
    path_point : ~pyarts3.arts.PropagationPathPoint, optional
        The propagation path point to use for the plot. Defaults to pos: [0, 0, 0], los: [0, 0].
    T_extrapolfac : float, optional
        Internal extrapolation factor
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig :
        As input if input.  Otherwise the created Figure.
    ax :
        As input if input.  Otherwise the created Axes.
    """

    if atm is None:
        ws = pyarts.Workspace()
        sp = []
        if spec1 != "AIR":
            sp.append(f"{spec1}")
        if spec2 != "AIR":
            sp.append(f"{spec2}")
        ws.abs_speciesSet(species=sp)
        basename = "planets/Earth/afgl/tropical/"
        toa = 1 + path_point.pos[0]
        ws.atm_fieldRead(toa=toa, basename=basename, missing_is_zero=1)
        atm = ws.atm_field(*path_point.pos)
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
                cia.spectral_propmat(
                    f=freqs, atm=atm, T_extrapolfac=T_extrapolfac, spec1=spec1, spec2=spec2)[:, 0],
                **kwargs
            )
    else:
        select_flat_ax(ax, 0).plot(
            freqs,
            data.spectral_propmat(
                f=freqs, atm=atm, T_extrapolfac=T_extrapolfac, spec1=spec1, spec2=spec2)[:, 0],
            **kwargs
        )

    return fig, ax
