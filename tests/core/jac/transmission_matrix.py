import pyarts3 as pyarts
import numpy as np


"""
Tests that the two-level transmission matrix has decent derivatives
for both distance and propagation matrix.  The test is done by perturbation.
"""

ArrayOfPropmatVector = pyarts.arts.ArrayOfPropmatVector
ArrayOfPropmatMatrix = pyarts.arts.ArrayOfPropmatMatrix
Propmat = pyarts.arts.Propmat
PropmatVector = pyarts.arts.PropmatVector
PropmatMatrix = pyarts.arts.PropmatMatrix
Vector = pyarts.arts.Vector
Tensor3 = pyarts.arts.Tensor3


def K(x, y, *, deriv=False):
    global v
    v[0] += 1
    if not deriv:
        return PropmatVector([v[0] * x + v[0] * 2 * y])
    else:
        return PropmatMatrix([[v[0]], [v[0] * 2], [0]])


def r(x, y, z, *, deriv=False):
    if not deriv:
        return z
    else:
        return [0, 0, 1]


x = np.array([1e-3, 1e-4, 1e-5])
y = np.array([1 / 30_000, 4e-5, 1e-2])
z = np.array([1000, 100, 10])

v = [0]
Ks = [K(x[i], y[i]) for i in range(len(x))]
v = [0]
dKs = [K(x[i], y[i], deriv=True) for i in range(len(x))]
rs = [r(x[i], y[i], z[i]) for i in range(len(x))]
drs = np.array([r(x[i], y[i], z[i], deriv=True) for i in range(len(x))])

Ks = ArrayOfPropmatVector(Ks)
dKs = ArrayOfPropmatMatrix(dKs)
rs = Vector(rs)
drs = Tensor3([drs * 0, drs])

T, dT = pyarts.arts.rtepack.two_level_exp(Ks, dKs, rs, drs)

dT = np.array(dT)

DERIV = 1e-8
dT2 = np.zeros_like(dT)
for i in range(1, len(x)):
    x2 = 1.0 * x
    y2 = 1.0 * y
    z2 = 1.0 * z
    x2[i - 1] += DERIV
    v = [0]
    Ks = [K(x2[i], y2[i]) for i in range(len(x))]
    v = [0]
    dKs = [K(x2[i], y2[i], deriv=True) for i in range(len(x))]
    rs = [r(x2[i], y2[i], z2[i]) for i in range(len(x))]
    drs = [r(x2[i], y2[i], z2[i], deriv=True) for i in range(len(x))]

    Ks = ArrayOfPropmatVector(Ks)
    dKs = ArrayOfPropmatMatrix(dKs)
    rs = Vector(rs)
    drs = Tensor3([drs, drs])

    T2, _ = pyarts.arts.rtepack.two_level_exp(Ks, dKs, rs, drs)

    dT2[i - 1, 0, 0] = ((np.array(T2) - np.array(T)) / DERIV)[i]

    x2 = 1.0 * x
    y2 = 1.0 * y
    z2 = 1.0 * z
    x2[i] += DERIV
    v = [0]
    Ks = [K(x2[i], y2[i]) for i in range(len(x))]
    v = [0]
    dKs = [K(x2[i], y2[i], deriv=True) for i in range(len(x))]
    rs = [r(x2[i], y2[i], z2[i]) for i in range(len(x))]
    drs = [r(x2[i], y2[i], z2[i], deriv=True) for i in range(len(x))]

    Ks = ArrayOfPropmatVector(Ks)
    dKs = ArrayOfPropmatMatrix(dKs)
    rs = Vector(rs)
    drs = Tensor3([drs, drs])

    T2, _ = pyarts.arts.rtepack.two_level_exp(Ks, dKs, rs, drs)

    dT2[i, 1, 0] = ((np.array(T2) - np.array(T)) / DERIV)[i]

    x2 = 1.0 * x
    y2 = 1.0 * y
    z2 = 1.0 * z
    y2[i - 1] += DERIV
    v = [0]
    Ks = [K(x2[i], y2[i]) for i in range(len(x))]
    v = [0]
    dKs = [K(x2[i], y2[i], deriv=True) for i in range(len(x))]
    rs = [r(x2[i], y2[i], z2[i]) for i in range(len(x))]
    drs = [r(x2[i], y2[i], z2[i], deriv=True) for i in range(len(x))]

    Ks = ArrayOfPropmatVector(Ks)
    dKs = ArrayOfPropmatMatrix(dKs)
    rs = Vector(rs)
    drs = Tensor3([drs, drs])

    T2, _ = pyarts.arts.rtepack.two_level_exp(Ks, dKs, rs, drs)

    dT2[i - 1, 0, 1] = ((np.array(T2) - np.array(T)) / DERIV)[i]

    x2 = 1.0 * x
    y2 = 1.0 * y
    z2 = 1.0 * z
    y2[i] += DERIV
    v = [0]
    Ks = [K(x2[i], y2[i]) for i in range(len(x))]
    v = [0]
    dKs = [K(x2[i], y2[i], deriv=True) for i in range(len(x))]
    rs = [r(x2[i], y2[i], z2[i]) for i in range(len(x))]
    drs = [r(x2[i], y2[i], z2[i], deriv=True) for i in range(len(x))]

    Ks = ArrayOfPropmatVector(Ks)
    dKs = ArrayOfPropmatMatrix(dKs)
    rs = Vector(rs)
    drs = Tensor3([drs, drs])

    T2, _ = pyarts.arts.rtepack.two_level_exp(Ks, dKs, rs, drs)

    dT2[i, 1, 1] = ((np.array(T2) - np.array(T)) / DERIV)[i]

    x2 = 1.0 * x
    y2 = 1.0 * y
    z2 = 1.0 * z
    z2[i - 1] += DERIV
    v = [0]
    Ks = [K(x2[i], y2[i]) for i in range(len(x))]
    v = [0]
    dKs = [K(x2[i], y2[i], deriv=True) for i in range(len(x))]
    rs = [r(x2[i], y2[i], z2[i]) for i in range(len(x))]
    drs = [r(x2[i], y2[i], z2[i], deriv=True) for i in range(len(x))]

    Ks = ArrayOfPropmatVector(Ks)
    dKs = ArrayOfPropmatMatrix(dKs)
    rs = Vector(rs)
    drs = Tensor3([drs, drs])

    T2, _ = pyarts.arts.rtepack.two_level_exp(Ks, dKs, rs, drs)

    dT2[i - 1, 0, 2] = ((np.array(T2) - np.array(T)) / DERIV)[i]

    x2 = 1.0 * x
    y2 = 1.0 * y
    z2 = 1.0 * z
    z2[i] += DERIV
    v = [0]
    Ks = [K(x2[i], y2[i]) for i in range(len(x))]
    v = [0]
    dKs = [K(x2[i], y2[i], deriv=True) for i in range(len(x))]
    rs = [r(x2[i], y2[i], z2[i]) for i in range(len(x))]
    drs = [r(x2[i], y2[i], z2[i], deriv=True) for i in range(len(x))]

    Ks = ArrayOfPropmatVector(Ks)
    dKs = ArrayOfPropmatMatrix(dKs)
    rs = Vector(rs)
    drs = Tensor3([drs, drs])

    T2, _ = pyarts.arts.rtepack.two_level_exp(Ks, dKs, rs, drs)

    dT2[i, 1, 2] = ((np.array(T2) - np.array(T)) / DERIV)[i]

assert np.allclose(dT2[np.where(dT2)] / dT[np.where(dT)], 1)
