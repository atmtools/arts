import matplotlib.pyplot as plt
import numpy as np


__all__ = ["default_fig_ax", "select_flat_ax", "stream_radiance"]


def default_fig_ax(fig=None, ax=None, nrows=1, ncols=1, N=-1, fig_kwargs={}, ax_kwargs={}):
    """Utility to create default matplotlib figure and axes if not provided.

    Parameters
    ----------
    fig : Figure, optional
        The matplotlib figure to draw on. Defaults to None for new figure.
    ax : Axes, optional
        The matplotlib axes to draw on. Defaults to None for new axes.
    fig_kwargs : dict, optional
        Keyword arguments for creating new figure if fig is None. Defaults to {}.
    ax_kwargs : dict, optional
        Keyword arguments for creating new axes if ax is None. Defaults to {}.

    Returns
    -------
    fig : Figure
        The matplotlib figure.
    ax : Axes
        The matplotlib axes.
    """
    fig = plt.figure(**fig_kwargs) if fig is None else fig
    if ax is None:
        ax = fig.subplots(nrows, ncols, **ax_kwargs)
        n = nrows * ncols
        if n > 1 and N < n and N >= 0:
            for a in ax.flatten()[N:]:
                fig.delaxes(a)

    return fig, ax


def select_flat_ax(ax, index):
    """Utility to select a single Axes from possibly multi-dimensional axes array.

    Parameters
    ----------
    ax : Axes or array of Axes
        The matplotlib axes or array of axes.
    index : int
        The index of the desired Axes in flattened order.

    Returns
    -------
    selected_ax : Axes
        The selected matplotlib axes.
    """
    if isinstance(ax, plt.Axes):
        return ax
    elif isinstance(ax, np.ndarray):
        return ax.flat[index]
    elif isinstance(ax, (list, tuple)):
        return ax[index]
    else:
        return ax.flatten()[index]


_STOKES = ["I", "Q", "U", "V"]


def stream_radiance(mu, up, down, *, fig=None, ax=None, stokes=None, **kwargs):
    """Draw radiance on the streams of a plane-parallel solver against the line-of-sight zenith angle.

    One panel per Stokes component.  Below 90 degrees the line of sight
    looks up and sees the downward radiance, above 90 degrees it looks down
    and sees the upward radiance.  Plotting several solvers on the same
    ``fig`` and ``ax`` compares them.

    Parameters
    ----------
    mu : array_like
        ``[stream]`` cosines of the streams in (0, 1], the same in both
        hemispheres, in any order.
    up : array_like
        ``[stream, stokes]`` radiance propagating upward.
    down : array_like
        ``[stream, stokes]`` radiance propagating downward.
    fig : Figure, optional
        The matplotlib figure to draw on.  Defaults to None for new figure.
    ax : array of Axes, optional
        One Axes per Stokes component to draw on.  Defaults to None for new
        axes.
    stokes : list of int, optional
        The Stokes components to show.  Defaults to all.
    **kwargs
        Passed to :func:`matplotlib.axes.Axes.plot`.

    Returns
    -------
    fig : Figure
        The matplotlib figure.
    ax : array of Axes
        The matplotlib axes, one per Stokes component.
    """
    mu = np.asarray(mu, dtype=float)
    up = np.asarray(up)
    down = np.asarray(down)
    stokes = list(range(up.shape[1])) if stokes is None else list(stokes)

    if fig is None:
        fig = plt.figure(figsize=(3.6 * len(stokes), 3.6), constrained_layout=True)
    if ax is None:
        ax = np.asarray(fig.subplots(1, len(stokes), sharex=True, squeeze=False))[0]
    ax = np.atleast_1d(ax)

    order = np.argsort(mu)
    mu = mu[order]
    # A NaN at 90 degrees separates the hemispheres, which are not continuous there
    za = np.concatenate((np.degrees(np.arccos(mu[::-1])), [90.0], 180.0 - np.degrees(np.arccos(mu))))
    for a, s in zip(ax, stokes):
        a.plot(za, np.concatenate((down[order, s][::-1], [np.nan], up[order, s])), **kwargs)
        a.set_title(f"Stokes {_STOKES[s]}")
        a.set_xlabel("Line-of-sight zenith angle [deg]")
        a.axvline(90.0, color="0.7", lw=0.5)
        a.grid(True, alpha=0.3)
    ax[0].set_ylabel("Radiance [W m$^{-2}$ Hz$^{-1}$ sr$^{-1}$]")
    if "label" in kwargs:
        ax[0].legend()
    return fig, ax
