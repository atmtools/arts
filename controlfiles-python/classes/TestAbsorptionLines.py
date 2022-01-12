import os
from pyarts.workspace import Workspace
from pyarts.classes import from_workspace
from pyarts.classes.AbsorptionLines import AbsorptionLines, ArrayOfAbsorptionLines
from pyarts.classes.quantum import QuantumIdentifier


ws = Workspace()
datapath = "../../../arts-cat-data/" if not os.getenv("ARTS_XML_DATA_DIR") else os.getenv("ARTS_XML_DATA_DIR")

# Init
al = AbsorptionLines()
aal = ArrayOfAbsorptionLines()
aaal = from_workspace(ws.abs_lines_per_species)

# Init is as expected
assert al.selfbroadening == False, "Bad init"
assert al.bathbroadening == False, "Bad init"
assert al.cutoff == "None", "Bad init"
assert al.mirroring == "None", "Bad init"
assert al.population == "LTE", "Bad init"
assert al.normalization == "None", "Bad init"
assert al.lineshapetype == "VP", "Bad init"
assert al.t0 == 296, "Bad init"
assert al.cutofffreq == -1, "Bad init"
assert al.linemixinglimit == -1, "Bad init"
assert al.quantumidentity == QuantumIdentifier(), "Bad init"
assert not al.broadeningspecies, "Bad init"
assert not al.lines, "Bad init"
assert al.OK, "Bad init"

# Read same file twice, in ARTS and external
ws.abs_speciesSet(species = ["O2-66"])
aal.readxml(os.path.join(datapath, "O2-66.xml"))
ws.abs_lines_per_speciesReadSpeciesSplitCatalog(basename = datapath)

# Everything should be the same (with silly not-empty test)
assert aal == aaal[0], "Bad load"
assert aal[0].lines[0].f0 != 0, "Bad frequency"
assert 2*aal[0].lines[0].f0 == 2*aaal[0][0].lines[0].f0, "Bad frequency"

al.set(aal[0])
assert al == aal[0], "Cannot set"

al2 = AbsorptionLines()
al.savexml("tmp.al.xml", "binary")
al2.readxml("tmp.al.xml")
assert al == al2
