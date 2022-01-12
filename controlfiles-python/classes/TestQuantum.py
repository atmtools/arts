from pyarts.workspace import Workspace
from pyarts.classes.quantum import QuantumIdentifier, ArrayOfQuantumIdentifier
from pyarts.classes.SpeciesTag import SpeciesTag
from pyarts.classes import from_workspace


# Get a workspace
ws = Workspace()

aqi = from_workspace(ws.band_identifiers)
qi = QuantumIdentifier()

qi.isot = "H2O-161"
qi.val["J"] = "J 3 2"

aqi.data = [qi]
assert qi == aqi[0]

aqi.size = 2
aqi[1].set(aqi[0])
assert aqi[0] == aqi[1]
aqi.append(qi)
assert aqi[0] == aqi[2]
for x in aqi:
    x.val["J"] = "J 3 2"

aqi2 = ArrayOfQuantumIdentifier()
aqi.savexml("tmp.aqi.xml", "ascii")
aqi2.readxml("tmp.aqi.xml")
assert aqi == aqi2
