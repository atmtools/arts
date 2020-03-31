from pyarts.workspace import Workspace
from pyarts.classes.EnergyLevelMap import EnergyLevelMap
from pyarts.classes import from_workspace


# Get a workspace
ws = Workspace()
datapath = "../"


elm1 = EnergyLevelMap()
elm1.readxml(datapath + "controlfiles/artscomponents/nlte/testdata/nlte_testdata.xml")
ws.ReadXML(ws.nlte_field , datapath + "controlfiles/artscomponents/nlte/testdata/nlte_testdata.xml")
elm2 = from_workspace(ws.nlte_field)

assert elm1, "Bad read"
assert elm1 == elm2, "Bad read"
assert elm1.data

elm3 = EnergyLevelMap()
elm3.set(elm1)

assert elm3 == elm2
