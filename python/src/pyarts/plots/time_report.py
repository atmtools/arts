import pyarts
import numpy as np
import matplotlib.pyplot as plt


def time_report(mode="f", *, clear=True, scale=1.0, fig=None):
    r = pyarts.arts.globals.time_report(clear)

    if fig is None:
        fig = plt.figure()

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

    if mode == "f":
        mapp = {}
        for x in r:
            for f in r[x]:
                if f not in mapp:
                    mapp[f] = []
                for p in r[x][f]:
                    tv = np.array([p[0].sec, p[1].sec])-m
                    tv *= scale
                    mapp[f].append((tv, x))
        
        keys = list(mapp.keys())
        keys.sort()
        
        ax = fig.add_subplot()
        
        for i in range(len(keys)):
            key = keys[i]
            for item in mapp[key]:
                tv, x = item
                ax.plot(tv, np.array([x, x]) / X + i - 0.5, "r")
        ax.set_yticks(range(len(keys)), keys)
    else:
        ax = fig.add_subplot()

        for x in r:
            for f in r[x]:
                for p in r[x][f]:
                    tv = np.array([p[0].sec, p[1].sec])-m
                    tv *= scale
                    ax.plot(tv, [x, x], "r", lw=3)

    return fig, ax, r
