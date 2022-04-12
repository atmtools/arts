import numpy as np

from pyarts.workspace import Workspace
from pyarts.classes.GriddedField1 import GriddedField1
from pyarts.classes.GriddedField2 import GriddedField2
from pyarts.classes.GriddedField3 import GriddedField3
from pyarts.classes.GriddedField4 import GriddedField4
from pyarts.classes.GriddedField5 import GriddedField5
from pyarts.classes.GriddedField6 import GriddedField6
from pyarts.classes.BasicTypes import String, ArrayOfString
from pyarts.classes.Vector import Vector
from pyarts.classes.Matrix import Matrix
from pyarts.classes.Tensor3 import Tensor3
from pyarts.classes.Tensor4 import Tensor4
from pyarts.classes.Tensor5 import Tensor5
from pyarts.classes.Tensor6 import Tensor6
from pyarts.classes import from_workspace


# Get a workspace
ws = Workspace()

# Get some variables
gf1 = from_workspace(ws.sideband_response)
gf2 = from_workspace(ws.surface_type_mask)
gf3 = from_workspace(ws.z_field_raw)
gf4 = from_workspace(ws.antenna_response)
gf5 = GriddedField5()
gf6 = GriddedField6()

assert not gf1
assert not gf2
assert not gf3
assert not gf4
assert not gf5
assert not gf6

assert gf1.OK
assert gf2.OK
assert gf3.OK
assert gf4.OK
assert gf5.OK
assert gf6.OK

gf1.grids = [Vector([1,2,3])]
assert not gf1.OK
gf1.data = Vector([5,4,3])
assert gf1.OK

gf2.grids = [Vector([1,2,3]), ArrayOfString([String("Hej"), String("hopp")])]
assert not gf2.OK
gf2.data = Matrix([[1, 2], [3, 4], [5, 6]])
assert gf2.OK

gf3.data = Tensor3(np.zeros((5, 3, 2)))
assert not gf3
assert not gf3.OK
gf3.grids = [Vector(np.linspace(0, 1, 5)), Vector(np.linspace(0, 1, 3)), Vector(Vector(np.linspace(0, 1, 2)))]
assert gf3.OK

gf4.data = Tensor4(np.zeros((5, 3, 2, 4)))
assert not gf4
assert not gf4.OK
gf4.grids = [Vector(np.linspace(0, 1, 5)), Vector(np.linspace(0, 1, 3)), Vector(np.linspace(0, 1, 2)), Vector(np.linspace(0, 1, 4))]
assert gf4.OK

gf5.data = Tensor5(np.zeros((5, 3, 2, 4)))
assert not gf5
assert not gf5.OK
gf5.grids = [Vector([0]), Vector(np.linspace(0, 1, 5)), Vector(np.linspace(0, 1, 3)), Vector(np.linspace(0, 1, 2)), Vector(np.linspace(0, 1, 4))]
assert gf5.OK

gf6.data = Tensor6(np.zeros((1, 2, 3, 4, 5, 6)))
assert not gf6
assert not gf6.OK
gf6.grids = [Vector([0]), Vector(np.linspace(0, 1, 2)), Vector(np.linspace(0, 1, 3)), Vector(np.linspace(0, 1, 4)), Vector(np.linspace(0, 1, 5)), Vector(np.linspace(0, 1, 6))]
assert gf6.OK

gf12 = GriddedField1()
gf22 = GriddedField2()
gf32 = GriddedField3()
gf42 = GriddedField4()
gf52 = GriddedField5()
gf62 = GriddedField6()

gf12.set(gf1)
gf22.set(gf2)
gf32.set(gf3)
gf42.set(gf4)
gf52.set(gf5)
gf62.set(gf6)

assert gf1 == gf12
assert gf2 == gf22
assert gf3 == gf32
assert gf4 == gf42
assert gf5 == gf52
assert gf6 == gf62

gf13 = GriddedField1()
gf23 = GriddedField2()
gf33 = GriddedField3()
gf43 = GriddedField4()
gf53 = GriddedField5()
gf63 = GriddedField6()

gf1.savexml("tmp.gf1.xml", "binary")
gf2.savexml("tmp.gf2.xml", "binary")
gf3.savexml("tmp.gf3.xml", "binary")
gf4.savexml("tmp.gf4.xml", "binary")
gf5.savexml("tmp.gf5.xml", "binary")
gf6.savexml("tmp.gf6.xml", "binary")

gf13.readxml("tmp.gf1.xml")
gf23.readxml("tmp.gf2.xml")
gf33.readxml("tmp.gf3.xml")
gf43.readxml("tmp.gf4.xml")
gf53.readxml("tmp.gf5.xml")
gf63.readxml("tmp.gf6.xml")

assert gf1 == gf13
assert gf2 == gf23
assert gf3 == gf33
assert gf4 == gf43
assert gf5 == gf53
assert gf6 == gf63
