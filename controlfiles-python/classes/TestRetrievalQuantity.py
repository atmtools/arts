from pyarts.workspace import Workspace
from pyarts.classes.QuantumIdentifier import QuantumIdentifier
from pyarts.classes.RetrievalQuantity import RetrievalQuantity, ArrayOfRetrievalQuantity
from pyarts.classes.SpeciesTag import SpeciesTag
from pyarts.classes.Vector import ArrayOfVector
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
rq.maintag = ""
rq.subtag = "HSE on"
rq.subsubtag = ""
rq.mode = ""
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

qn = QuantumIdentifier("Transition")
aost = [SpeciesTag(), SpeciesTag()]
aost[0].setFromString("H2O-161")
aost[1].setFromString("H2O-MPM89")

rq2 = RetrievalQuantity()
rq2.maintag = ""
rq2.subtag = ""
rq2.subsubtag = ""
rq2.mode = ""
rq2.analytical = 1
rq2.target.type = "Special"
rq2.target.subtype = "SurfaceString"
rq2.target.perturbation = 0.0
rq2.grids = ArrayOfVector([p, lat, lon])
rq2.target.string_key = "I am a string"
rq2.target.quantumidentity = qn
rq2.target.specieslist = aost

arq.append(rq2)

rq2.target.type = "Line"
rq2.target.subtype = "VMR"

arq.append(rq2)

rq2.target.type = "Special"
rq2.target.subtype = "ArrayOfSpeciesTagVMR"

arq.append(rq2)

arq2 = ArrayOfRetrievalQuantity()
arq.savexml("tmp.arq.xml", "ascii")
arq2.readxml("tmp.arq.xml")

try:
    assert arq == arq2
except:
    ws.ReadXML(ws.jacobian_quantities, "tmp.arq.xml")
    assert arq == arq2
    Warning("We had a failure that should not be!!!")
