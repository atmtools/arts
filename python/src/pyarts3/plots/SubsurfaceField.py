import numpy as np
import matplotlib.pyplot as plt

def plot(
    subsurf_field,
    *,
    fig=None,
    alts=np.linspace(0, 1e5),
    lats=0,
    lons=0,
    ygrid=None,
    keys=None,
):
    """Plot the subsurface field parameters in a default manner.

    Args:
        subsurf_field (pyarts3.arts.SubsurfaceField): A subsurface field
        fig (optional): The matplotlib figure to draw on. Defaults to None for new figure.
        alts (optional): A grid to plot on - must after broadcast with lats and lons be 1D. Defaults to np.linspace(0, 1e5).
        lats (optional): A grid to plot on - must after broadcast with alts and lons be 1D. Defaults to 0.
        lons (optional): A grid to plot on - must after broadcast with alts and lats be 1D. Defaults to 0.
        ygrid (optional): Choice of y-grid for plotting.  Uses broadcasted alts if None. Defaults to None.
        keys (optional): A list of keys to plot. Defaults to None for all keys in keys().

    Returns:
        fig: as input or a new figure
        subs: list of subplots
    """
    alts, lats, lons = np.broadcast_arrays(alts, lats, lons)
    v = subsurf_field(alts, lats, lons)

    keys = v[0].keys() if keys is None else keys
    N = len(keys)
    n = int(np.ceil(np.sqrt(N))) + 1

    if fig is None:
        fig = plt.figure(figsize=(5 * n, 5 * n))

    subs = []
    for i in range(N):
        subs.append(fig.add_subplot(n, n, i + 1))
        subs[-1].plot(
            [x[keys[i]] for x in v],
            alts if ygrid is None else ygrid,
            label=keys[i],
        )
        subs[-1].legend()
    return fig, subs
