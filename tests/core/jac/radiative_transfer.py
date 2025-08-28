import pyarts3 as pyarts
import numpy as np


"""
Tests that the two-level radiative transfer produces decent
results.  The test is done by perturbation. 
"""


ArrayOfPropmatVector = pyarts.arts.ArrayOfPropmatVector
ArrayOfPropmatMatrix = pyarts.arts.ArrayOfPropmatMatrix
ArrayOfStokvecVector = pyarts.arts.ArrayOfStokvecVector
ArrayOfStokvecMatrix = pyarts.arts.ArrayOfStokvecMatrix
Propmat = pyarts.arts.Propmat
PropmatVector = pyarts.arts.PropmatVector
PropmatMatrix = pyarts.arts.PropmatMatrix
StokvecVector = pyarts.arts.StokvecVector
StokvecMatrix = pyarts.arts.StokvecMatrix
Vector = pyarts.arts.Vector
Tensor3 = pyarts.arts.Tensor3


def K(x, *, deriv=False):
    global v
    v[0] += 1
    if not deriv:
        return PropmatVector([v[0] * x])
    else:
        return PropmatMatrix([[v[0]], [0], [0]])


def r(z, *, deriv=False):
    global v
    v[0] += 1
    if not deriv:
        return v[0] * z
    else:
        return [0, 0, v[0]]


def J(y, deriv=False):
    global v
    v[0] += 1
    if not deriv:
        return StokvecVector([v[0] * y])
    else:
        return StokvecMatrix([[0], [v[0]], [0]])


x = np.array([1e-2, 1e-3, 1e-4, 5e-4, 1e-4])
y = np.array([3000, 400, 100, 300, 200])
z = np.array([1000, 100, 1000 / 3, 33, 5])

N = len(x)

v = [0]
Ks = [K(x[i]) for i in range(len(x))]
v = [0]
dKs = [K(x[i], deriv=True) for i in range(len(x))]
v = [0]
rs = [r(z[i]) for i in range(len(x))]
v = [0]
drs = np.array([r(z[i], deriv=True) for i in range(len(x))])
v = [0]
Js = ArrayOfStokvecVector([J(y[i]) for i in range(len(x))])
v = [0]
dJs = ArrayOfStokvecMatrix([J(y[i], deriv=True) for i in range(len(x))])

Ks = ArrayOfPropmatVector(Ks)
dKs = ArrayOfPropmatMatrix(dKs)
rs = Vector(rs)
drs = Tensor3([drs * 0, drs])

T, dT = pyarts.arts.rtepack.two_level_exp(Ks, dKs, rs, drs)

I0 = StokvecVector([300])
I, dI = pyarts.arts.rtepack.two_level_radiative_transfer(T, dT, Js, dJs, I0)

dI = np.array(dI)

DERIVS = [x[0] / 1e8, y[0] / 1e8, z[0] / 1e8]
dI2 = np.zeros((N, 3, 4))
for j in range(len(x)):
    x2 = 1.0 * x
    y2 = 1.0 * y
    z2 = 1.0 * z
    x2[j] += DERIVS[0]

    v = [0]
    Ks = [K(x2[i]) for i in range(len(x))]
    v = [0]
    dKs = [K(x2[i], deriv=True) for i in range(len(x))]
    v = [0]
    rs = [r(z2[i]) for i in range(len(x))]
    v = [0]
    drs = [r(z2[i], deriv=True) for i in range(len(x))]
    v = [0]
    Js = [J(y2[i]) for i in range(len(x))]
    v = [0]
    dJs = [J(y2[i], deriv=True) for i in range(len(x))]

    Ks = ArrayOfPropmatVector(Ks)
    dKs = ArrayOfPropmatMatrix(dKs)
    rs = Vector(rs)
    drs = Tensor3([drs, drs])

    T2, dT2 = pyarts.arts.rtepack.two_level_exp(Ks, dKs, rs, drs)
    I2, _ = pyarts.arts.rtepack.two_level_radiative_transfer(
        T2, dT2, Js, dJs, I0
    )

    dI2[j, 0] = (I2 - I) / DERIVS[0]

    x2 = 1.0 * x
    y2 = 1.0 * y
    z2 = 1.0 * z
    y2[j] += DERIVS[1]

    v = [0]
    Ks = [K(x2[i]) for i in range(len(x))]
    v = [0]
    dKs = [K(x2[i], deriv=True) for i in range(len(x))]
    v = [0]
    rs = [r(z2[i]) for i in range(len(x))]
    v = [0]
    drs = [r(z2[i], deriv=True) for i in range(len(x))]
    v = [0]
    Js = [J(y2[i]) for i in range(len(x))]
    v = [0]
    dJs = [J(y2[i], deriv=True) for i in range(len(x))]

    Ks = ArrayOfPropmatVector(Ks)
    dKs = ArrayOfPropmatMatrix(dKs)
    rs = Vector(rs)
    drs = Tensor3([drs, drs])

    T2, dT2 = pyarts.arts.rtepack.two_level_exp(Ks, dKs, rs, drs)
    I2, _ = pyarts.arts.rtepack.two_level_radiative_transfer(
        T2, dT2, Js, dJs, I0
    )

    dI2[j, 1] = (I2 - I) / DERIVS[1]

    x2 = 1.0 * x
    y2 = 1.0 * y
    z2 = 1.0 * z
    z2[j] += DERIVS[2]

    v = [0]
    Ks = [K(x2[i]) for i in range(len(x))]
    v = [0]
    dKs = [K(x2[i], deriv=True) for i in range(len(x))]
    v = [0]
    rs = [r(z2[i]) for i in range(len(x))]
    v = [0]
    drs = [r(z2[i], deriv=True) for i in range(len(x))]
    v = [0]
    Js = [J(y2[i]) for i in range(len(x))]
    v = [0]
    dJs = [J(y2[i], deriv=True) for i in range(len(x))]

    Ks = ArrayOfPropmatVector(Ks)
    dKs = ArrayOfPropmatMatrix(dKs)
    rs = Vector(rs)
    drs = Tensor3([drs, drs])

    T2, dT2 = pyarts.arts.rtepack.two_level_exp(Ks, dKs, rs, drs)
    I2, _ = pyarts.arts.rtepack.two_level_radiative_transfer(
        T2, dT2, Js, dJs, I0
    )

    dI2[j, 2] = (I2 - I) / DERIVS[2]

assert np.allclose(dI2[np.where(dI2)] / dI[np.where(dI)], 1)
