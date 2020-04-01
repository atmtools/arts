from pyarts.workspace import Workspace
from pyarts.classes.SpeciesTag import SpeciesTag
from pyarts.classes import from_workspace


# Get a workspace
ws = Workspace()

aast = from_workspace(ws.abs_species)
st = SpeciesTag("H2O-161")
ws.abs_speciesSet(species=["H2O-161"])

assert st == aast[0][0]
