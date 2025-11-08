import pyarts3 as pyarts

import numpy as np

vec = np.random.random(size=3)

def saveit():
  ws = pyarts.Workspace(False)

  ws.atm_fieldInit(toa=100e3, default_isotopologue="None")
  ws.atm_fieldIGRF()

  #ws.savexml("field_with_op.xml")

  #return ws.atm_field["mag_u"].data(*vec)

def readit():
  pass
  #ws = pyarts.Workspace.fromxml("field_with_op.xml")

  #return ws.atm_field["mag_u"].data(*vec)

#v1 = saveit()
#v2 = readit()

#assert v1 == v2

#print(v1, 'and', v2, 'are equal!')
