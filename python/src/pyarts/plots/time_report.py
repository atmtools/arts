import pyarts
import numpy as np
import matplotlib.pyplot as plt


def time_report(mode="f", *, clear=True, scale=1.0):
    r = pyarts.arts.globals.time_report(clear)

    m = None
    for x in r:
        for f in r[x]:
            for p in r[x][f]:
                if m is None:
                    m = p[0].sec
                m = min(m, p[0].sec)
                m = min(m, p[1].sec)

    if mode == "f":
        X = 0
        mapp = {}
        for x in r:
            X = max(x, X)
            for f in r[x]:
                if f not in mapp:
                    mapp[f] = []
                for p in r[x][f]:
                    tv = np.array([p[0].sec, p[1].sec])-m
                    tv *= scale
                    mapp[f].append((tv, x))
        
        keys = list(mapp.keys())
        keys.sort()
        for i in range(len(keys)):
            key = keys[i]
            for item in mapp[key]:
                tv, x = item
                plt.plot(tv, np.array([x, x]) / X + i - 0.5, "r")
        plt.yticks(range(len(keys)), keys)
    else:
        for x in r:
            for f in r[x]:
                for p in r[x][f]:
                    tv = np.array([p[0].sec, p[1].sec])-m
                    tv *= scale
                    plt.plot(tv, [x, x], "r", lw=3)
