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

# Empty NumPy arrays must retain their dynamic dimensions when converted to
# rtepack arrays, while their trailing dimensions form one fixed-size block.
rtepack_empty_arrays = [
    (pyarts.arts.StokvecMatrix, (0, 2, 4), float),
    (pyarts.arts.StokvecTensor3, (0, 2, 3, 4), float),
    (pyarts.arts.PropmatMatrix, (0, 2, 7), float),
    (pyarts.arts.MuelmatMatrix, (0, 2, 4, 4), float),
    (pyarts.arts.MuelmatTensor3, (0, 2, 3, 4, 4), float),
    (pyarts.arts.MuelmatTensor4, (0, 2, 3, 4, 4, 4), float),
    (pyarts.arts.MuelmatTensor5, (0, 2, 3, 4, 5, 4, 4), float),
    (pyarts.arts.SpecmatTensor3, (0, 2, 3, 4, 4), complex),
]

for rtepack_type, shape, dtype in rtepack_empty_arrays:
    converted = rtepack_type(np.empty(shape, dtype=dtype))
    result = np.asarray(converted)
    assert result.shape == shape
    assert result.dtype == np.dtype(dtype)

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
