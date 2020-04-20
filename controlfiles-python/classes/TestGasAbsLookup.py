from pyarts.workspace import Workspace
from pyarts.classes.GasAbsLookup import GasAbsLookup
from pyarts.classes import from_workspace


# Get a workspace
ws = Workspace()


gal = from_workspace(ws.abs_lookup)
assert isinstance(gal, GasAbsLookup)

gal2 = GasAbsLookup()
gal.savexml("tmp.gal.xml", "binary")
gal2.readxml("tmp.gal.xml")
assert gal == gal2
