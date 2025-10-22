""" Plotting routine for the sensor response """

import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt

__all__ = [
    'plot',
]


def plot(measurement_sensor: pyarts.arts.ArrayOfSensorObsel, *, fig=None, keys: str | list = "f", pol: str | pyarts.arts.Stokvec = "I", **kwargs):
    """Plot the sensor observational element array.

    .. note::
        Any option other than "f" will sort the sensor observational element array by the
        corresponding key. The object lives via a pointer, so the original
        measurement_sensor is modified in-place.

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
                                            frequency_grid = np.linspace(-50e6, 50e6, 101))
        pyarts.plots.ArrayOfSensorObsel.plot(ws.measurement_sensor)

    Parameters
    ----------
    measurement_sensor : ~pyarts3.arts.ArrayOfSensorObsel
        A sensor observation element array.
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    subs : Subaxis, optional
        List of subplots to add to. Defaults to None for a new subplot.
    keys : str | list
        The keys to use for plotting. Options are in :class:`~pyarts3.arts.SensorKeyType`.
    pol : str | pyarts3.arts.Stokvec
        The polarization to use for plotting. Defaults to "I", constructs a :class:`~pyarts3.arts.Stokvec`.

    Returns
    -------
    fig : As input
        As input.
    subs : As input
        As input.
    """

    if isinstance(keys, str):
        keys = [keys]
    N = len(keys)

    if N == 0:
        return fig, None

    if fig is None:
        fig = plt.figure(figsize=(10 * N, 10))

    pol = pyarts.arts.Stokvec(pol)

    map = {
        pyarts.arts.SensorKeyType.f: None,
        pyarts.arts.SensorKeyType.alt: 0,
        pyarts.arts.SensorKeyType.lat: 1,
        pyarts.arts.SensorKeyType.lon: 2,
        pyarts.arts.SensorKeyType.za: 3,
        pyarts.arts.SensorKeyType.aa: 4,
    }

    subs = []
    for isub in range(N):
        key = pyarts.arts.SensorKeyType(keys[isub])

        i = map[key]

        subs.append(fig.add_subplot(1, N, isub + 1))
        for elem in measurement_sensor:
            if i is None:
                v = elem.weight_matrix.reduce(pol, along_poslos=True)
            else:
                v = elem.weight_matrix.reduce(pol, along_freq=True)

            x = elem.f_grid if i is None else elem.poslos[:, i]

            if len(x) == 1:
                subs[-1].plot(x, v, marker="o", linestyle="None")
            else:
                subs[-1].plot(x, v)

    return fig, subs
