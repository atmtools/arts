import numpy as np
import pyarts3 as pyarts


def not_owns_data(nd):
    assert not nd.flags["OWNDATA"]


def owns_data(nd):
    assert nd.flags["OWNDATA"]


def writable(nd):
    assert nd.flags["WRITEABLE"]


def not_writable(nd):
    assert not nd.flags["WRITEABLE"]


ts = [pyarts.arts.Vector,
      pyarts.arts.Matrix,
      pyarts.arts.Tensor3,
      pyarts.arts.Tensor4,
      pyarts.arts.Tensor5,
      pyarts.arts.Tensor6,
      pyarts.arts.Tensor7]

for i in range(7):
    v = ts[i](np.linspace(0, 1, 2**(i+1)).reshape(*([2] * (i + 1))))

    a = np.array(v)
    b = np.array(v, copy=False)
    c = np.array(v, copy=False, dtype=float)
    d = np.array(v, copy=True)
    e = np.array(v, copy=True, dtype=float)

    try:
        f = np.array(v, copy=False, dtype=int)
        exit(1)
    except ValueError:
        pass

    test1 = [owns_data(x) for x in [a, d, e]]
    test2 = [not_owns_data(x) for x in [b, c]]
    test3 = [writable(x) for x in [a, b, c, d, e]]

v = pyarts.arts.StokvecMatrix([[1, 2, 3], [4, 5, 6]])

a = np.array(v)
b = np.array(v, copy=False)
c = np.array(v, copy=False, dtype=float)
d = np.array(v, copy=True)
e = np.array(v, copy=True, dtype=float)

try:
    f = np.array(v, copy=False, dtype=int)
    exit(1)
except ValueError:
    pass

test1 = [owns_data(x) for x in [a, d, e]]
test2 = [not_owns_data(x) for x in [b, c]]
test3 = [writable(x) for x in [a, b, c, d, e]]

v = pyarts.arts.DescendingGrid([3, 2, 1])

a = np.array(v)
b = np.array(v, copy=False)
c = np.array(v, copy=False, dtype=float)
d = np.array(v, copy=True)
e = np.array(v, copy=True, dtype=float)

try:
    f = np.array(v, copy=False, dtype=int)
    exit(1)
except ValueError:
    pass

test1 = [owns_data(x) for x in [a, d, e]]
test2 = [not_owns_data(x) for x in [b, c]]
test3 = [writable(x) for x in [a, d, e]]
test4 = [not_writable(x) for x in [b, c]]

ts = [pyarts.arts.AscendingGrid,
      pyarts.arts.LatGrid,
      pyarts.arts.LonGrid,
      pyarts.arts.ZenGrid,
      pyarts.arts.AziGrid]

for t in ts:
    v = t([1, 2, 3])

    a = np.array(v)
    b = np.array(v, copy=False)
    c = np.array(v, copy=False, dtype=float)
    d = np.array(v, copy=True)
    e = np.array(v, copy=True, dtype=float)

    try:
        f = np.array(v, copy=False, dtype=int)
        exit(1)
    except ValueError:
        pass

    test1 = [owns_data(x) for x in [a, d, e]]
    test2 = [not_owns_data(x) for x in [b, c]]
    test3 = [writable(x) for x in [a, d, e]]
    test4 = [not_writable(x) for x in [b, c]]
