import numpy as np
import scipy as sp
from pyarts.arts import Propmat as Propmat
from pyarts.arts import PropmatVector as PropmatVector
from pyarts.arts import PropmatMatrix as PropmatMatrix
from pyarts.arts import ArrayOfPropmatVector as ArrayOfPropmatVector
from pyarts.arts import ArrayOfPropmatMatrix as ArrayOfPropmatMatrix
from pyarts.arts import Vector as Vector
from pyarts.arts import Tensor3 as Tensor3
from pyarts.arts.rtepack import two_level_exp


a, b, c, d, u, v, w = np.random.random((7))
print(f"a, b, c, d, u, v, w = {a, b, c, d, u, v, w}")
a += 1


def Ks(x):
    offset = 0

    K = Propmat([offset + a * x, b * x, c * x, d * x, u * x, v * x, w * x])
    dK = Propmat([a, b, c, d, u, v, w])

    delta = 1e-8
    x += delta
    K2 = Propmat([offset + a * x, b * x, c * x, d * x, u * x, v * x, w * x])

    return K, K2, dK, (K2 - K) / delta, delta


def mat(dK):
    x = np.empty((1, 1, 7))

    x[0, 0] = dK

    return PropmatMatrix(x)


xs = np.logspace(-3, -1, 51)
t1 = []
t2 = []
fail = []
for x in xs:
    K, K2, dK, dK2, delta = Ks(x)
    
    T = sp.linalg.expm(-K.as_matrix())
    
    KS = ArrayOfPropmatVector([PropmatVector([K]), PropmatVector([K])])
    KS2 = ArrayOfPropmatVector([PropmatVector([K2]), PropmatVector([K2])])
    DKS = ArrayOfPropmatMatrix([mat(dK), mat(dK)])
    RS = Vector([1, 1])
    DRS = Tensor3(np.zeros((2, 2, 1)))
    
    TS, DTS = two_level_exp(KS, DKS, RS, DRS)
    TS2, DTS2 = two_level_exp(KS2, DKS, RS, DRS)
    
    assert np.allclose((TS2[1] - TS[1]) * 1e8 / DTS[0][0] / 2, 1, rtol=1e-3)
