from pyarts.workspace import Workspace
from pyarts.classes.RetrievalQuantity import RetrievalQuantity, ArrayOfRetrievalQuantity
from pyarts.classes.Vector import Vector, ArrayOfVector
from pyarts.classes import from_workspace


# Get a workspace
ws = Workspace()
ws.jacobianInit()
ws.atmosphere_dim = 3
p = from_workspace(ws.p_grid)
lat = from_workspace(ws.lat_grid)
lon = from_workspace(ws.lon_grid)
p.data = [3, 2, 1]
lat.data = [1, 2, 3, 4]
lon.data = [1, 2, 3, 4, 5]
ws.jacobianAddTemperature(g1=p.data, g2=lat.data, g3=lon.data)

arq = from_workspace(ws.jacobian_quantities)
rq = RetrievalQuantity()
rq.maintag = "Atmospheric temperatures"
rq.subtag = "HSE on"
rq.subsubtag = "From propagation matrix"
rq.mode = "abs"
rq.analytical = 1
rq.target.type = "Atm"
rq.target.subtype = "Temperature"
rq.target.perturbation = 0.1
rq.grids = ArrayOfVector([p, lat, lon])

assert rq == arq[0]
arq.size = 2
arq[1].set(rq)
assert arq[1] == arq[0]
arq.append(rq)
assert arq[1] == arq[2]

arq2 = ArrayOfRetrievalQuantity()
arq.savexml("tmp.arq.xml", "ascii")
arq2.readxml("tmp.arq.xml")

try:
    assert arq == arq2
except:
    ws.ReadXML(ws.jacobian_quantities, "tmp.arq.xml")
    assert arq == arq2
    Warning("We had a failure that should not be!!!")
