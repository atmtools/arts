from pyarts.workspace import Workspace
from pyarts.classes.PropagationMatrix import PropagationMatrix, ArrayOfPropagationMatrix
from pyarts.classes.StokesVector import StokesVector, ArrayOfStokesVector
from pyarts.classes import from_workspace


# Get a workspace
ws = Workspace()

# Get the agenda
apm = from_workspace(ws.dpropmat_clearsky_dx)
asv = from_workspace(ws.dnlte_source_dx)
pm = PropagationMatrix(1.0, 4, 4, 4, 4)
sv = StokesVector(1.0, 4, 4, 4, 4)

assert isinstance(apm, ArrayOfPropagationMatrix)
apm.size = 6
apm[0].setData(1, 4, 3, 2, 0.0)
apm[1].setData(1, 4, 2, 4, 1.0)
apm[2].setData(2, 3, 1, 8, 2.0)
apm[3].setData(3, 2, 2, 16, 3.14)
apm[4].setData(5, 1, 1, 32, 6.28)
apm[4].set(apm[3])
apm[5].set(pm)

assert not apm[0]
assert apm[1]
assert apm[3] == apm[4]
assert apm[5] == pm

assert isinstance(asv, ArrayOfStokesVector)
asv.size = 6
asv[0].setData(1, 4, 3, 2, 0.0)
asv[1].setData(1, 4, 2, 4, 1.0)
asv[2].setData(2, 3, 1, 8, 2.0)
asv[3].setData(3, 2, 2, 16, 3.14)
asv[4].setData(5, 1, 1, 32, 6.28)
asv[5].set(sv)
asv[4].set(asv[3])
assert not asv[0]
assert asv[1]
assert asv[5] == sv

apm.savexml("tmp.apm.xml", "binary")
asv.savexml("tmp.asv.xml", "binary")
apm2 = ArrayOfPropagationMatrix()
asv2 = ArrayOfStokesVector()
apm2.readxml("tmp.apm.xml")
asv2.readxml("tmp.asv.xml")
assert apm2 == apm
assert asv2 == asv
