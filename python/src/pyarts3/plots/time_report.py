import pyarts3 as pyarts
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


def time_report(*, mode="plot", clear=True, scale=1.0, fig=None, mintime=None):
    """Plots the time report.

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
        The figure to draw on in a plotting mode.  By default creates a new figure.

    Return
    ------
    r : dict
        The time report
    *out : tuple
        If mode is "plot", returns the matplotlib figure and axis.
        If mode is "table", returns a string containing the time report in Markdown table format.
    """

    r = pyarts.arts.globals.time_report(clear)

    m = None
    X = 0
    for x in r:
        X = max(x, X)
        for f in r[x]:
            for p in r[x][f]:
                if m is None:
                    m = p[0].timestamp()
                m = min(m, p[0].timestamp())
                m = min(m, p[1].timestamp())

    dt = {}
    res = {}
    for x in r:
        for f in r[x]:
            if f not in res:
                res[f] = []
                dt[f] = [1e309, 0.0]
            for p in r[x][f]:
                tv = np.array([p[0].timestamp(), p[1].timestamp()]) - m
                tv *= scale
                res[f].append([tv, np.array([x, x])])
                dt[f][0] = np.min([tv[0], dt[f][0]])
                dt[f][1] += tv[1] - tv[0]

    if mintime is not None:
        keys = np.array(list(res.keys()))
        for key in keys:
            if dt[key][1] < mintime:
                del res[key]
                del dt[key]

    unit = time_report_unit(scale)

    if mode == "plot":
        out = time_report_plot(res, dt, unit, fig)
    elif mode == "table":
        out = [time_report_table(res, dt, unit)]
    else:
        raise ValueError(
            f"Unknown mode {mode}, see method description for valid modes."
        )

    return r, *out


def time_report_unit(scale):
    if scale == 1.0:
        return "s"
    elif scale == 1e3:
        return "ms"
    elif scale == 1e6:
        return "µs"
    else:
        return f"{scale} x s"


def time_report_plot(res, dt, unit, fig):
    """Plots the time report.

    Parameters
    ----------
    res : dict
        The flattened time report.
    dt : dict
        The time report with the minimum time for each method.
    unit : str
        The unit of the time axis, e.g., 's', 'ms', 'µs', etc.
    fig : matplotlib figure, optional
        The figure to draw on.

    Returns
    -------
    fig : matplotlib figure
    ax : matplotlib axis on figure
    """

    keys = np.array(list(res.keys()))
    vs = np.array([dt[key][0] for key in keys])
    keys = keys[np.argsort(vs)]

    if fig is None:
        fig = plt.figure()

    ax = fig.add_subplot()

    colors = cm.get_cmap("viridis", len(keys))

    i = 0
    for key in keys:
        v = res[key]
        plt.plot(v[0][0], v[0][1], color=colors.colors[i], lw=3, label=key)

        if len(v) > 0:
            for x in v[1:]:
                plt.plot(*x, color=colors.colors[i], lw=3)
        i += 1

    ax.legend(ncols=4, loc="center left", bbox_to_anchor=(1, 0.5))
    ax.set_ylabel("Core ID")
    ax.set_xlabel(f"Time [{unit}]")

    return fig, ax


def time_report_table(res, dt, unit):
    """Prints the time report in text format.

    Parameters
    ----------
    res : dict
        The flattened time report.
    dt : dict
        The time report with the minimum time for each method.

    Returns
    -------
    str
        A string containing the time report in Markdown table format.
    """

    keys = np.array(list(res.keys()))
    vs = np.array([dt[key][1] for key in keys])
    keys = keys[np.argsort(vs)][::-1]

    out = f"| Method    | Total Time [{unit}]  | Min Time [{unit}] | Max Time [{unit}] | Average Time [{unit}] | Times Called |\n"
    out += "| --------- | ------------ | ------------ | -------- | ------------ | ------------ |\n"

    for key in keys:
        numc = len(res[key])
        avgt = dt[key][1] / numc
        dts = [v[0][1] - v[0][0] for v in res[key]]
        mint = min(dts)
        maxt = max(dts)
        tott = dt[key][1]
        out += f"| {key} | {round(tott,1)} | {round(mint,1)} | {round(maxt,1)} | {round(avgt,1)} | {numc} |\n"

    return out
