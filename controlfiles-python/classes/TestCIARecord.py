import os

from pyarts.workspace import Workspace
from pyarts.classes.CIARecord import CIARecord, ArrayOfCIARecord
from pyarts.classes import from_workspace


# Get a workspace
ws = Workspace()
datapath = "../../arts-xml-data/" if not os.getenv("ARTS_XML_DATA_DIR") else os.getenv("ARTS_XML_DATA_DIR")
fn = os.path.join(datapath, "spectroscopy/cia/borysow/Borysow_CIA_JUICE_SWI.xml")

acr1 = ArrayOfCIARecord()
assert not acr1, "Bad init"

acr1.readxml(fn)
ws.ReadXML(ws.abs_cia_data , fn)
acr2 = from_workspace(ws.abs_cia_data)

assert acr1, "Bad read"
assert acr1 == acr2, "Bad read"
assert isinstance(acr1[0], CIARecord), "Bad type"

acr3 = ArrayOfCIARecord()
acr3.set(acr1)
assert acr1 == acr3, "Bad read"
