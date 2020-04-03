import os

from pyarts.workspace import Workspace
from pyarts.classes.ScatteringMetaData import ScatteringMetaData
from pyarts.classes.SingleScatteringData import SingleScatteringData
from pyarts.classes import from_workspace


# Get a workspace
ws = Workspace()
datapath = "../../arts-xml-data/" if not os.getenv("ARTS_XML_DATA_DIR") else os.getenv("ARTS_XML_DATA_DIR")
fn1 = os.path.join(datapath, 'scattering/H2O_ice/MieSphere_R1.00000e+00um.meta.xml')
fn2 = os.path.join(datapath, 'scattering/H2O_ice/MieSphere_R1.00000e+00um.xml')

smd1 = ScatteringMetaData()
smd1.readxml(fn1)
smd2 = from_workspace(ws.scat_meta_single)
ws.ReadXML(ws.scat_meta_single, fn1)
smd3 = ScatteringMetaData()
smd3.set(smd2)

assert smd1 == smd2
assert smd1 == smd3

ssd1 = SingleScatteringData()
ssd1.readxml(fn2)
ssd2 = from_workspace(ws.scat_data_single)
ws.ReadXML(ws.scat_data_single, fn2)
ssd3 = SingleScatteringData()
ssd3.set(ssd2)

assert ssd1 == ssd2
assert ssd1 == ssd3

smd4 = ScatteringMetaData()
smd1.savexml("tmp.smd.xml", "binary")
smd4.readxml("tmp.smd.xml")
assert smd1 == smd4

ssd4 = SingleScatteringData()
ssd1.savexml("tmp.ssd.xml", "binary")
ssd4.readxml("tmp.ssd.xml")
assert ssd1 == ssd4
