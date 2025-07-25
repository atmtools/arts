import pyarts
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


def time_report(*, clear=True, scale=1.0, fig=None):
    """ Plots the time report.

    The time report is available only when
    ARTS has been compiled with the CMake option ``-DENABLE_ARTS_PROFILING=OFF``,
    which is not part of the default distribution.

    Nevertheless, there is enough helpful parts in being able to see the
    state of parallelism and report the time of internal methods that
    this method is part of the distribution.

    Parameters
    ----------
    clear : bool, optional
        Whether or not to clear the time-report.  Default is True.
    scale : float, optional
        The scale of the time axis, defaults to 1.0.
    fig : matplotlib figure, optional
        The figure to draw on.  By default creates a new figure.
    
    Return
    ------
    fig : matplotlib figure
    ax : matplotlib axis on figure
    r : dict
        The time report

    """
    r = pyarts.arts.globals.time_report(clear)

    if fig is None:
        fig = plt.figure()

    ax = fig.add_subplot()

    m = None
    X = 0
    for x in r:
        X = max(x, X)
        for f in r[x]:
            for p in r[x][f]:
                if m is None:
                    m = p[0].sec
                m = min(m, p[0].sec)
                m = min(m, p[1].sec)

    dt = {}
    res = {}
    for x in r:
        for f in r[x]:
            if f not in res:
                res[f] = []
                dt[f] = 1e309
            for p in r[x][f]:
                tv = np.array([p[0].sec, p[1].sec])-m
                tv *= scale
                res[f].append([tv, np.array([x, x])])
                dt[f] = np.min([tv[0], dt[f]])

    keys = np.array(list(res.keys()))
    vs = np.array([dt[key] for key in keys])
    keys = keys[np.argsort(vs)]


    colors = cm.get_cmap('viridis', len(keys))

    i = 0
    for key in keys:
        v = res[key]
        plt.plot(v[0][0], v[0][1], color=colors.colors[i], lw=3, label=key)

        if len(v) > 0:
            for x in v[1:]:
                plt.plot(*x, color=colors.colors[i], lw=3)
        i += 1

    ax.legend(ncols= 4, loc='center left', bbox_to_anchor=(1, 0.5))
    ax.set_ylabel("Core ID")
    if scale == 1.0:
        ax.set_xlabel("Time [s]")
    elif scale == 1e3:
        ax.set_xlabel("Time [ms]")
    elif scale == 1e6:
        ax.set_xlabel("Time [µs]")
    else:
        ax.set_xlabel(f"Time [{scale} x s]")

    return fig, ax, r
