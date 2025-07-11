import pyarts

import numpy as np

vec = np.random.random(size=3)

def saveit():
  ws = pyarts.Workspace(False)

  ws.atmospheric_fieldInit(toa=100e3, default_isotopologue="None")
  ws.atmospheric_fieldIGRF()

  ws.savexml("field_with_op.xml")

  return ws.atmospheric_field["mag_u"].data(*vec)

def readit():
  ws = pyarts.Workspace.fromxml("field_with_op.xml")

  return ws.atmospheric_field["mag_u"].data(*vec)

v1 = saveit()
v2 = readit()

assert v1 == v2

print(v1, 'and', v2, 'are equal!')
