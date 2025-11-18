import pyarts3 as pyarts
import numpy as np

x = pyarts.arts.ArrayOfIndex()

assert np.allclose(x, [])

x.append(0)

assert np.allclose(x, [0]), "Cannot append"

pyarts.arts.Index(1).savexml("index1.xml")

x.appendxml("index1.xml")

assert np.allclose(x, [0] + [1]), "Cannot appendxml"

pyarts.arts.ArrayOfIndex([2, 3]).savexml("index2.xml")

x.extendxml("index2.xml")

assert np.allclose(x, [0] + [1] + [2, 3]), "Cannot extendxml"

x.extend([4, 5, 6])

assert np.allclose(x, [0] + [1] + [2, 3] + [4, 5, 6]), "Cannot extend"
