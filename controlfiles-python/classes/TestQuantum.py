from pyarts.workspace import Workspace
from pyarts.classes.QuantumIdentifier import QuantumIdentifier, ArrayOfQuantumIdentifier
from pyarts.classes.SpeciesTag import SpeciesTag
from pyarts.classes import from_workspace


# Get a workspace
ws = Workspace()

aqi = from_workspace(ws.band_identifiers)
qi = QuantumIdentifier()

qi.spec_ind = 0
qi.type = "None"
qi.type = "All"
qi.type = "EnergyLevel"
qi.type = "Transition"
qi.upp["J"] = 3
qi.low["J"] = 2

aqi.data = [qi]
assert qi == aqi[0]

aqi.size = 2
aqi[1].set(aqi[0])
assert aqi[0] == aqi[1]
aqi.append(qi)
assert aqi[0] == aqi[2]
for x in aqi:
    x.upp["J"] = 3
    x.low["J"] = 2

aqi2 = ArrayOfQuantumIdentifier()
aqi.savexml("tmp.aqi.xml", "ascii")
aqi2.readxml("tmp.aqi.xml")
assert aqi == aqi2
