import pyarts
import numpy as np
import matplotlib.pyplot as plt


def plot(measurement_sensor, *, fig=None, keys="f", pol="I"):
    """Plot the sensor obsel array in a default manner.

    .. note::
        Any option other than "f" will sort the sensor obsel array by the
        corresponding key. The object lives via a pointer, so the original
        measurement_sensor is modified in-place.

        This should not massively affect any radiative transfer computations,
        but it might result in different results down on the floating point
        epsilon error level.

    Args:
        measurement_sensor (pyarts.arts.MeasurementSensor): A sensor observation element array
        fig (optional): The matplotlib figure to draw on. Defaults to None for new figure.
        subs (optional): List of subplots to add to. Defaults to None for a new subplot.
        keys (str, list): The keys to use for plotting. Options are in pyarts.arts.SensorKeyType.
        pol (str): The polarization to use for plotting. Defaults to "I", constructs a pyarts.arts.Stokvec.
    Returns:
        fig: as input or a new figure
        subs: list of subplots
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
            v = np.einsum(
                "ijk,k->j" if i is None else "ijk,k->i", elem.weight_matrix, pol
            )
            x = elem.f_grid if i is None else elem.poslos[:, i]

            if len(x) == 1:
                subs[-1].plot(x, v, marker="o", linestyle="None")
            else:
                subs[-1].plot(x, v)

    return fig, subs
