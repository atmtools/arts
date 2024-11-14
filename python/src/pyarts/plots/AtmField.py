import numpy as np
import matplotlib.pyplot as plt

def plot(
    atm_field,
    *,
    fig=None,
    alts=np.linspace(0, 1e5),
    lats=0,
    lons=0,
    ygrid=None,
    keep_basic=True,
    keep_specs=True,
    keep_isots=False,
    keep_nlte=False,
    keep_ssprops=True,
):
    """Plot the atmospheric field parameters in a default manner.

    Args:
        atm_field (pyarts.arts.AtmField): An atmospheric field
        fig (optional): The matplotlib figure to draw on. Defaults to None for new figure.
        alts (optional): A grid to plot on - must after broadcast with lats and lons be 1D. Defaults to np.linspace(0, 1e5).
        lats (optional): A grid to plot on - must after broadcast with alts and lons be 1D. Defaults to 0.
        lons (optional): A grid to plot on - must after broadcast with alts and lats be 1D. Defaults to 0.
        ygrid (optional): Choice of y-grid for plotting.  Uses broadcasted alts if None. Defaults to None.
        keep_basic (bool, optional): Forwarded to pyarts.arts.AtmPoint::keys. Defaults to True.
        keep_specs (bool, optional): Forwarded to pyarts.arts.AtmPoint::keys. Defaults to True.
        keep_isots (bool, optional): Forwarded to pyarts.arts.AtmPoint::keys. Defaults to False.
        keep_nlte (bool, optional): Forwarded to pyarts.arts.AtmPoint::keys. Defaults to False.
        keep_ssprops (bool, optional): Forwarded to pyarts.arts.AtmPoint::keys. Defaults to True.

    Returns:
        fig: as input or a new figure
        subs: list of subplots
    """
    alts, lats, lons = np.broadcast_arrays(alts, lats, lons)
    v = atm_field.at(alts, lats, lons)

    keys = v[0].keys(
        keep_basic=keep_basic,
        keep_specs=keep_specs,
        keep_isots=keep_isots,
        keep_nlte=keep_nlte,
        keep_ssprops=keep_ssprops,
    )
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
