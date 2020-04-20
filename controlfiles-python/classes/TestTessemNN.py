from pyarts.workspace import Workspace
from pyarts.classes.TessemNN import TessemNN
from pyarts.classes import from_workspace


# Get a workspace
ws = Workspace()

tnn = from_workspace(ws.tessem_netv)
tnn.nb_inputs = 1
tnn.nb_outputs = 2
tnn.nb_cache = 3
tnn.b1 = 4
tnn.b2 = 5
tnn.w1 = 6
tnn.w2 = 7
tnn.x_min = 8
tnn.x_max = 9
tnn.y_min = 10
tnn.y_max = 11

assert tnn == tnn
tnn2 = from_workspace(ws.tessem_neth)
tnn2.set(tnn)
assert tnn == tnn2

# tnn3 = TessemNN()
# tnn.savexml("tmp.tnn.xml", "binary")
# tnn3.readxml("tmp.tnn.xml")
# assert tnn == tnn3
