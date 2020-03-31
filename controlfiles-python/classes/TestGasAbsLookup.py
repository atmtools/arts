from pyarts.workspace import Workspace
from pyarts.classes.GasAbsLookup import GasAbsLookup
from pyarts.classes import from_workspace


# Get a workspace
ws = Workspace()
fn = "../controlfiles/testdata/testdoit_gas_abs_lookup.xml"

gal1 = GasAbsLookup()
assert not gal1, "Bad init"

gal1.readxml(fn)
ws.ReadXML(ws.abs_lookup , fn)
gal2 = from_workspace(ws.abs_lookup)

assert gal1, "Bad read"
assert gal1 == gal2, "Bad read"
