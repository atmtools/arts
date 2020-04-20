from pyarts.workspace import Workspace
from pyarts.classes.SpeciesTag import SpeciesTag, ArrayOfArrayOfSpeciesTag
from pyarts.classes import from_workspace


# Get a workspace
ws = Workspace()

aast = from_workspace(ws.abs_species)
st = SpeciesTag("H2O-161")
ws.abs_speciesSet(species=["H2O-161"])
aast.size = 2
aast[1].size = 1
aast[1][0].set(aast[0][0])

assert st == aast[0][0]
assert st == aast[1][0]

aast.savexml("tmp.aast.xml", "binary")
aast2 = ArrayOfArrayOfSpeciesTag()
aast2.readxml("tmp.aast.xml")
assert aast == aast2
