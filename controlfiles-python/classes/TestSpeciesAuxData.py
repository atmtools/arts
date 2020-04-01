import os

from pyarts.workspace import Workspace
from pyarts.classes.SpeciesAuxData import SpeciesAuxData
from pyarts.classes import from_workspace


# Get a workspace
ws = Workspace()
datapath = "../../arts-xml-data/" if not os.getenv("ARTS_XML_DATA_DIR") else os.getenv("ARTS_XML_DATA_DIR")
fn = os.path.join(datapath, 'planets/Mars/isotopratio_Mars.xml')

sad1 = SpeciesAuxData()
sad1.readxml(fn)
sad2 = from_workspace(ws.isotopologue_ratios)
ws.ReadXML(ws.isotopologue_ratios, fn)
sad3 = SpeciesAuxData()
sad3.set(sad2)

assert sad1 == sad2
assert sad1 == sad3
