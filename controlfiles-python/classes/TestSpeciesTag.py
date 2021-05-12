from pyarts.workspace import Workspace
from pyarts.classes.SpeciesTag import (SpeciesTag, ArrayOfSpeciesTag,
                                       ArrayOfArrayOfSpeciesTag)
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

# Test ArrayOfSpeciesTag
ws.ArrayOfSpeciesTagCreate("ast")

ws.ArrayOfSpeciesTagSet(ws.ast, "O2-66, O2-68")
assert str(ws.ast.value) == '"O2-66-*-*, O2-68-*-*"'

ws.ArrayOfSpeciesTagSet(ws.ast, ["O2-66", "O2-68"])
assert str(ws.ast.value) == '"O2-66-*-*, O2-68-*-*"'

ws.ArrayOfSpeciesTagSet(ws.ast, "O2-66")
assert str(ws.ast.value) == '"O2-66-*-*"'

ws.ArrayOfSpeciesTagSet(ws.ast, "")
assert str(ws.ast.value) == '""'

ws.ArrayOfSpeciesTagSet(ws.ast, None)
assert str(ws.ast.value) == '""'

ws.ArrayOfSpeciesTagSet(ws.ast, ArrayOfSpeciesTag([SpeciesTag("O2-66")]))
assert str(ws.ast.value) == '"O2-66-*-*"'
