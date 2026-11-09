""" Plotting routine for the sensor response """

import numpy
import matplotlib
import pyarts3 as pyarts
from .common import default_fig_ax, select_flat_ax

__all__ = [
    'plot',
]


def plot(data: pyarts.arts.ArrayOfSensor,
         *,
         fig: matplotlib.figure.Figure | None = None,
         ax: matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes] | None = None,
         keys: str | list = "f",
         pol: str | pyarts.arts.Stokvec = "I",
         **kwargs) -> tuple[matplotlib.figure.Figure, matplotlib.axes.Axes | list[matplotlib.axes.Axes] | numpy.ndarray[matplotlib.axes.Axes]]:
    """Plot the observation elements of all the sensors.

    .. rubric:: Example

    .. plot::
        :include-source:

        import pyarts3 as pyarts
        import numpy as np

        ws = pyarts.Workspace()
        ws.measurement_sensorsSimpleGaussian(std = 10e6, pos = [100e3, 0, 0], los = [180.0, 0.0],
                                            freq_grid = np.linspace(-50e6, 50e6, 101))
        pyarts.plots.ArrayOfSensor.plot(ws.measurement_sensors)

    Parameters
    ----------
    data : ~pyarts3.arts.ArrayOfSensor
        The sensors.
    fig : ~matplotlib.figure.Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : ~matplotlib.axes.Axes | list[~matplotlib.axes.Axes] | ~numpy.ndarray[~matplotlib.axes.Axes] | None, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    keys : str | list
        The keys to use for plotting. Options are in :class:`~pyarts3.arts.SensorKeyType`.
    pol : str | ~pyarts3.arts.Stokvec
        The polarization to use for plotting. Defaults to "I", constructs a :class:`~pyarts3.arts.Stokvec`.
    **kwargs : keyword arguments
        Additional keyword arguments to pass to the plotting functions.

    Returns
    -------
    fig : matplotlib.figure.Figure
        As input if input.  Otherwise the created Figure.
    ax : matplotlib.axes.Axes
        As input if input.  Otherwise the created Axes.
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
        pyarts.arts.SensorKeyType.zen: 3,
        pyarts.arts.SensorKeyType.azi: 4,
    }

    for isub in range(N):
        key = pyarts.arts.SensorKeyType(keys[isub])

        i = map[key]
        for sensor in data:
            for weights in sensor.elements:
                v = weights.reduce(pol,
                                   along_poslos=i is None,
                                   along_freq=i is not None)

                x = sensor.freq_grid if i is None else sensor.geometry[:, i]

                select_flat_ax(ax, isub).plot(x, v, **kwargs)

    return fig, ax
