import pyarts.pyarts_cpp as cxx

def oi(ws):
    print("oi")
    ws.atmosphere_dim = 3

x = cxx.CallbackFunction(oi)

ws = cxx.Workspace()
assert not ws.atmosphere_dim.init
x(ws)
assert ws.atmosphere_dim.init
assert ws.atmosphere_dim.value.val == 3
