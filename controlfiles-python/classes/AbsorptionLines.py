import os
from pyarts.workspace import Workspace
from pyarts.classes import from_workspace
from pyarts.classes.AbsorptionLines import AbsorptionLines, ArrayOfAbsorptionLines


ws = Workspace()
datapath = "../../arts-xml-data/" if not os.getenv("ARTS_XML_DATA_DIR") else os.getenv("ARTS_XML_DATA_DIR")
fn1 = os.path.join(datapath, 'spectroscopy/cat/')
fn2 = os.path.join(datapath, 'spectroscopy/cat/O2-66.xml')

# Init
al = AbsorptionLines()
aal = ArrayOfAbsorptionLines()
aaal = from_workspace(ws.abs_lines_per_species)

# Init is as expected
assert al.selfbroadening == False, "Bad init"
assert al.bathbroadening == False, "Bad init"
assert al.cutoff == 0, "Bad init"
assert al.mirroring == 0, "Bad init"
assert al.population == 0, "Bad init"
assert al.normalization == 0, "Bad init"
assert al.lineshapetype == 0, "Bad init"
assert al.t0 == 296, "Bad init"
assert al.cutofffreq == -1, "Bad init"
assert al.linemixinglimit == -1, "Bad init"
assert al.quantumidentity.spec == 0, "Bad init"
assert al.quantumidentity.isot == 0, "Bad init"
assert al.quantumidentity.type == 3, "Bad init"
assert len(al.quantumidentity.upperqn.data) == len(al.quantumidentity.lowerqn.data), "Bad init"
assert not any(al.quantumidentity.upperqn.data), "Bad init"
assert not any(al.quantumidentity.lowerqn.data), "Bad init"
assert not al.localquantumnumbers, "Bad init"
assert not al.broadeningspecies, "Bad init"
assert not al.lines, "Bad init"
assert al.OK, "Bad init"

# Read same file twice, in ARTS and external
ws.abs_speciesSet(species = ["O2-66"])
aal.readxml(fn2)
ws.abs_lines_per_speciesReadSpeciesSplitCatalog(basename = fn1)

# Everythin should be the same (with silly not-empty test)
assert aal == aaal[0], "Bad load"
assert aal[0].lines[0].f0 != 0, "Bad frequency"
assert 2*aal[0].lines[0].f0 == 2*aaal[0][0].lines[0].f0, "Bad frequency"

del al
del aal
del aaal
del ws
