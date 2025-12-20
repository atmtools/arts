import numpy as np
import pyarts3 as pyarts
import scipy as sp


for i in range(1000):
    K = pyarts.arts.Propmat(np.random.normal(0.0, 10.0, size=7))

    rK = sp.linalg.sqrtm(K.as_matrix())
    rKp = pyarts.arts.rtepack.sqrt(K)

    assert np.allclose(rKp, rK)
    assert np.allclose(rKp @ rKp, K.as_matrix())

for i in range(7):
    K = pyarts.arts.Propmat([0, 0, 0, 0, 0, 0, 0])
    K[i] = 2

    rK = sp.linalg.sqrtm(K.as_matrix())
    rKp = pyarts.arts.rtepack.sqrt(K)

    assert np.allclose(rKp, rK)
    assert np.allclose(rKp @ rKp, K.as_matrix())
