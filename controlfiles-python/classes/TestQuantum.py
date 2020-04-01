from pyarts.workspace import Workspace
from pyarts.classes.QuantumIdentifier import QuantumIdentifier
from pyarts.classes.SpeciesTag import SpeciesTag
from pyarts.classes import from_workspace


# Get a workspace
ws = Workspace()

aqi = from_workspace(ws.band_identifiers)
qi = QuantumIdentifier()

qi.spec = 0
qi.isot = 0
qi.type = 3
qi.type = 2
qi.type = 1
qi.type = 0
qi.upperqn["J"] = 3
qi.lowerqn["J"] = 2

aqi.data = [qi]
assert qi == aqi[0]

aqi.size = 2
aqi[1].set(aqi[0])
assert aqi[0] == aqi[1]
aqi.append(qi)
assert aqi[0] == aqi[2]
