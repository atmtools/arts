import pyarts.pyarts_cpp as cxx
import pyarts

def oi(ws):
    print("oi")
    ws.atmosphere_dim = 3

x = cxx.CallbackFunction(oi)

ws = pyarts.workspace.Workspace()
assert not ws.atmosphere_dim.init
x(ws)
assert ws.atmosphere_dim.init
assert ws.atmosphere_dim.value.val == 3
