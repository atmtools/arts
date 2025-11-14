""" Plotting routine for the sensor response """

import pyarts3 as pyarts
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.ArrayOfSensorObsel,
         *,
         fig=None,
         ax=None,
         keys: str | list = "f",
         pol: str | pyarts.arts.Stokvec = "I",
         **kwargs):
    """Plot the sensor observational element array.

    .. note::
        Any option other than "f" will sort the sensor observational element array by the
        corresponding key. The object lives via a pointer, so the original
        data is modified in-place.

        This should not massively affect any radiative transfer computations,
        but it might result in different results down on the floating point
        epsilon error level.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        ws = pyarts.Workspace()
        ws.measurement_sensorSimpleGaussian(std = 10e6, pos = [100e3, 0, 0], los = [180.0, 0.0],
                                            freq_grid = np.linspace(-50e6, 50e6, 101))
        pyarts.plots.ArrayOfSensorObsel.plot(ws.measurement_sensor)

    Parameters
    ----------
    data : ~pyarts3.arts.ArrayOfSensorObsel
        A sensor observation element array.
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        Not used (function creates its own subplots). Accepted for API consistency.
    keys : str | list
        The keys to use for plotting. Options are in :class:`~pyarts3.arts.SensorKeyType`.
    pol : str | ~pyarts3.arts.Stokvec
        The polarization to use for plotting. Defaults to "I", constructs a :class:`~pyarts3.arts.Stokvec`.
    **kwargs : keyword arguments
        Additional keyword arguments passed to the plotting function.

    Returns
    -------
    fig : As input
        As input.
    ax : list
        List of matplotlib axes objects.
    """

    keys = [keys] if isinstance(keys, str) else keys
    N = len(keys)

    if N == 0:
        return fig, ax

    fig, ax = default_fig_ax(fig, ax, 1, N, fig_kwargs={'figsize': (10 * N, 10)})

    pol = pyarts.arts.Stokvec(pol)

    map = {
        pyarts.arts.SensorKeyType.f: None,
        pyarts.arts.SensorKeyType.alt: 0,
        pyarts.arts.SensorKeyType.lat: 1,
        pyarts.arts.SensorKeyType.lon: 2,
        pyarts.arts.SensorKeyType.za: 3,
        pyarts.arts.SensorKeyType.aa: 4,
    }

    for isub in range(N):
        key = pyarts.arts.SensorKeyType(keys[isub])

        i = map[key]
        for elem in data:
            v = elem.weight_matrix.reduce(pol,
                                          along_poslos=i is None,
                                          along_freq=i is not None)

            x = elem.f_grid if i is None else elem.poslos[:, i]

            select_flat_ax(ax, isub).plot(x, v, **kwargs)

    return fig, ax
